"""Export sequences and family IDs only. Held-out assay responses never leave the core."""
import argparse
import hashlib
import json
import shutil
from common import ROOT, rows, settings, table, atomic_text, sha256, guard


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    cfg = settings()
    split = list(rows(ROOT / "config/split_manifest.tsv"))
    sources = {r["tool"]: r for r in rows(ROOT / "results/upstream_software.tsv")}
    ds = {r["target"]: r for r in rows(ROOT / "config/datasets.tsv")}
    records = []
    for method in cfg["models"]:
        src = sources[method]
        for target in sorted(ds):
            for strategy in ["family", "random"]:
                group = [r for r in split if r["target"] == target and r["strategy"] == strategy and float(r["threshold"]) == cfg["family_threshold"]]
                dest = ROOT / f"remote/bundles/{method}_{target}_{strategy}"
                dest.mkdir(parents=True, exist_ok=True)
                train = [dict(sequence=r["sequence"], family_id=r["family_id"]) for r in group if r["split"] == "train"]
                test = [dict(sequence=r["sequence"], family_id=r["family_id"]) for r in group if r["split"] == "test"]
                table(dest / "train.tsv", train)
                table(dest / "test_sequences.tsv", test)
                job = dict(method=method, target=target, strategy=strategy, threshold=cfg["family_threshold"],
                     source_repository=src["source"], source_commit=src["source_commit"], software_licence=src["licence"],
                     licence_status="approved" if src["licence"] != "unclear" else "blocked", training_seed=cfg["master_seed"], generation_seeds=cfg["generation_seeds"],
                     generation_budget=cfg["generation_budget"], maximum_draws=cfg["generation_max_draws"],
                     minimum_length=int(ds[target]["minimum_length"]), maximum_length=int(ds[target]["maximum_length"]), motif_length=int(ds[target]["nominal_length"]),
                     forbidden_primers=[ds[target]["fold_forward"], ds[target]["fold_reverse"]],
                     configuration_hash=hashlib.sha256(json.dumps(cfg,sort_keys=True).encode()).hexdigest(),
                     input_hashes={n:sha256(dest/n) for n in ["train.tsv", "test_sequences.tsv"]},
                     context="unsupervised_or_selex_only_generation", assay_values_used=False,
                     architecture="upstream_CNN_PHMM_VAE" if method == "raptgen" else "upstream_conditional_diffusion_plus_VAE",
                     epochs=1000, early_stopping_patience=50, internal_validation="10_percent_remaining_training_families_label_blind",
                     generation_settings="standard_normal_2D_prior; upstream_PHMM_sampler" if method == "raptgen" else "blocked_no_authorized_implementation", keep_checkpoints=False)
                atomic_text(dest / "job.json", json.dumps(job,indent=2,sort_keys=True))
                records.append(dict(method=method,target=target,strategy=strategy,source_commit=src["source_commit"],
                     bundle=str(dest.relative_to(ROOT)), job_sha256=sha256(dest/"job.json"),
                     training_data_sha256=sha256(dest/"train.tsv"),test_sequences_sha256=sha256(dest/"test_sequences.tsv"),
                     configuration_hash=job["configuration_hash"], status="prepared_not_executed" if method == "raptgen" else "licence_blocked"))
    table(ROOT / "results/gpu_jobs.tsv", records)
    print("Prepared hashed GPU inputs without assay labels or responses; no GPU model run claimed")


if __name__ == "__main__":
    main()
