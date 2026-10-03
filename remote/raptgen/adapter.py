"""Call the pinned RaptGen architecture, loss, trainer and PHMM sampler directly."""
import argparse
import csv
from collections import defaultdict
import gzip
import hashlib
import json
from pathlib import Path
import random
import subprocess
import time
from datetime import datetime, timezone


def digest(path):
    h = hashlib.sha256()
    with open(path,"rb") as f:
        for block in iter(lambda:f.read(1048576),b""):
            h.update(block)
    return h.hexdigest()


def read(path):
    with open(path) as f:
        return list(csv.DictReader(f,delimiter="\t"))


def write(path, records, fields):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path,"wt") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(records)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--bundle",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    a = p.parse_args()
    job = json.loads((a.bundle/"job.json").read_text())
    for name, expected in job["input_hashes"].items():
        if digest(a.bundle/name) != expected:
            raise ValueError("Training input hash mismatch")
    if job["assay_values_used"] or job["licence_status"] != "approved":
        raise ValueError("Unapproved job")
    import numpy as np
    import torch
    from raptgen.models import CNN_PHMM_VAE, train
    from raptgen.data import one_hot_index, ProfileHMMSampler
    if not torch.cuda.is_available():
        raise RuntimeError("GPU is required for this hosted job")
    if torch.cuda.get_device_capability()[0] > 7:
        raise RuntimeError("Published CUDA 10.2/PyTorch 1.5 stack requires compatible older GPU such as T4; no modern-GPU compatibility claimed")
    start, wall = datetime.now(timezone.utc).isoformat(), time.monotonic()
    torch.manual_seed(job["training_seed"])
    torch.cuda.manual_seed_all(job["training_seed"])
    np.random.seed(job["training_seed"])
    random.seed(job["training_seed"])
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False
    records = read(a.bundle/"train.tsv")
    fids = sorted({r["family_id"] for r in records},key=lambda f:hashlib.sha256((str(job["training_seed"])+f).encode()).digest())
    validation = set(fids[:max(1,len(fids)//10)])
    class LengthBatches:
        def __init__(self, recs):
            self.dataset = recs
        def __iter__(self):
            buckets = defaultdict(list)
            for r in self.dataset:
                buckets[len(r["sequence"])].append(one_hot_index(r["sequence"]))
            batches = []
            for values in buckets.values():
                random.shuffle(values)
                chunks = [values[i:i+256] for i in range(0,len(values),256)]
                if len(chunks)>1 and len(chunks[-1]) == 1:
                    chunks[-2].extend(chunks.pop())
                for chunk in chunks:
                    if len(chunk)<2:
                        # BatchNorm needs two examples. Repetition does not introduce a new sequence.
                        chunk = chunk*2
                    batches.append(torch.tensor(chunk,dtype=torch.long))
            random.shuffle(batches)
            yield from batches
    train_rows = [r for r in records if r["family_id"] not in validation]
    val_rows = [r for r in records if r["family_id"] in validation]
    if not train_rows or not val_rows:
        raise RuntimeError("Internal family validation split is empty")
    model = CNN_PHMM_VAE(motif_len=job["motif_length"],embed_size=2).cuda()
    checkpoint = a.output/"model.mdl"
    losses = train(epochs=job["epochs"], model=model,train_loader=LengthBatches(train_rows),test_loader=LengthBatches(val_rows),optimizer=torch.optim.Adam(model.parameters()),
                   device="cuda",save_dir=a.output,model_str="model.mdl",threshold=job["early_stopping_patience"],beta_schedule=True,force_matching=True,force_epochs=50,logs=False)
    if not checkpoint.exists():
        raise RuntimeError("Upstream training produced no finite checkpoint")
    model.load_state_dict(torch.load(checkpoint))
    model.eval()
    def score(s):
        with torch.no_grad():
            (t,e),_,_ = model(torch.tensor([one_hot_index(s)],dtype=torch.long,device="cuda"),deterministic=True)
        sampler = ProfileHMMSampler(t[0].cpu().numpy(),e[0].cpu().numpy(),proba_is_log=True)
        return float(sampler.calc_seq_proba(s))/len(s)
    write(a.output/"heldout_scores.tsv",[dict(sequence=r["sequence"],score=score(r["sequence"]),score_type="posterior_mean_conditional_PHMM_log_probability_per_nt") for r in read(a.bundle/"test_sequences.tsv")], ["sequence","score","score_type"])
    generation_stats = []
    for seed in job["generation_seeds"]:
        np.random.seed(seed)
        torch.manual_seed(seed)
        seen, generated, draws, valid = set(), [], 0, 0
        while len(generated)<job["generation_budget"] and draws<job["maximum_draws"]:
            with torch.no_grad():
                t,e = model.decoder(torch.randn(128,2,device="cuda"))
            for tr,em in zip(t.cpu().numpy(),e.cpu().numpy()):
                draws += 1
                sampler = ProfileHMMSampler(tr,em,proba_is_log=True)
                s = sampler.sample(sequence_only=True).replace("_","").upper()
                if not job["minimum_length"] <= len(s) <= job["maximum_length"] or set(s)-set("ACGT") or any(pr in s for pr in job["forbidden_primers"]):
                    continue
                valid += 1
                if s not in seen:
                    seen.add(s)
                    generated.append(dict(sequence=s))
                if len(generated)==job["generation_budget"]:
                    break
        if len(generated)<job["generation_budget"]:
            raise RuntimeError("Generation unique budget unattained; seed must be reported as failed")
        write(a.output/f"generations_{seed}.tsv.gz",generated,["sequence"])
        generation_stats.append(dict(seed=seed,raw=draws,valid=valid,unique=len(generated)))
    outputs = ["heldout_scores.tsv"]+[f"generations_{s}.tsv.gz" for s in job["generation_seeds"]]
    manifest = {k:job[k] for k in ["method","target","strategy","source_repository","source_commit","configuration_hash","input_hashes","context","training_seed","generation_seeds"]}
    manifest.update(dict(status="complete",assay_values_used=False,architecture=job["architecture"],generation_settings=job["generation_settings"],
        package_versions=subprocess.check_output(["python","-m","pip","freeze"],universal_newlines=True).splitlines(), cuda_version=torch.version.cuda,
        gpu_model=torch.cuda.get_device_name(),gpu_memory_bytes=torch.cuda.get_device_properties(0).total_memory,
        peak_allocated_gpu_memory_bytes=torch.cuda.max_memory_allocated(),start_time=start,end_time=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-wall,
        epochs_executed=len(losses),checkpoint_hash=digest(checkpoint),checkpoint_location="regenerate_from_job_and_source_commit; not_retained",generation_statistics=generation_stats,output_hashes={n:digest(a.output/n) for n in outputs}))
    (a.output/"run_manifest.json").write_text(json.dumps(manifest,indent=2))
    checkpoint.unlink()


if __name__ == "__main__":
    main()
