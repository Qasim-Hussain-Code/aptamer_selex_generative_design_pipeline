"""Download one round, validate, process, verify compact output, then remove raw input."""
import argparse
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
from common import ROOT, rows, table, dataset, download, guard, sha256, failure, disk_bytes, atomic_text


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", choices=["all", "tg2", "integrin"], default="all")
    a = p.parse_args()
    manifest = list(rows(ROOT / "results/selex_manifest.tsv"))
    inputs = [ROOT / "config/datasets.tsv", ROOT / "results/selex_manifest.tsv"]
    config_hash = sha256(inputs[0])
    selected = [r for r in manifest if a.target == "all" or r["target"] == a.target]
    pilot = None
    for r in selected:
        target, rnd = r["target"], r["round"]
        out = ROOT / f"data/processed/{target}/{rnd}.counts.tsv.gz"
        qc = ROOT / f"results/qc/{target}_{rnd}.tsv"
        state = ROOT / f"logs/state/fetch_{target}_{rnd}.json"
        if out.exists() and qc.exists() and state.exists():
            saved = json.loads(state.read_text())
            if saved == dict(config_hash=config_hash, provider_md5=r["md5"], output_hash=sha256(out), qc_hash=sha256(qc)):
                pilot = next(rows(qc))
                print(f"fetch {target} round {rnd}: counts and input metadata validated, skipping", flush=True)
                continue
        raw = ROOT / "data/raw" / Path(r["url"]).name
        extra = int(r["bytes"]) * 8
        guard(extra)
        if pilot:
            ratio = int(pilot["compact_bytes"]) / int(pilot["compressed_raw_bytes"])
            remaining = sum(int(x["bytes"]) for x in selected)
            projected_extra = remaining * ratio + max(int(x["bytes"]) for x in selected) * (1 + int(pilot["sqlite_peak_bytes"])/int(pilot["compressed_raw_bytes"]))
            projected = guard(projected_extra)
            table(ROOT / "results/disk_projection.tsv", [dict(pilot_target=pilot["target"], pilot_round=pilot["round"],
                  raw_bytes=pilot["compressed_raw_bytes"], compact_bytes=pilot["compact_bytes"], compact_to_raw_ratio=ratio,
                  projected_peak_bytes=projected, safety_margin_fraction=0.1, policy="one_raw_file_at_a_time; no_checkpoints")])
        try:
            download(r["url"], raw, r["run"], r["md5"], r["bytes"])
            subprocess.run([sys.executable, str(ROOT / "scripts/resource_wrapper.py"), "--stage", f"preprocess_{target}_{rnd}", "--", sys.executable,
                            str(ROOT / "scripts/05_preprocess_selex.py"), "--input", str(raw), "--target", target, "--round", rnd], check=True)
            pilot = next(rows(qc))
            if int(pilot["raw_reads"]) != int(r["expected_reads"]):
                raise ValueError("Read count differs from provider metadata")
            atomic_text(state, json.dumps(dict(config_hash=config_hash, provider_md5=r["md5"], output_hash=sha256(out), qc_hash=sha256(qc))))
        except Exception as e:
            failure("fetch", r["run"], e, "delete raw/partial; preserve validated compact outputs", target)
            raise
        finally:
            raw.unlink(missing_ok=True)
            raw.with_name(raw.name + ".part").unlink(missing_ok=True)
    table(ROOT / "results/preprocessing_counts.tsv", [next(rows(p)) for p in sorted((ROOT / "results/qc").glob("*.tsv"))])


if __name__ == "__main__":
    main()
