"""Reject mismatched provenance, path escapes, altered hashes and incomplete seed budgets."""
import argparse
import json
import re
import shutil
from pathlib import Path
from common import ROOT, rows, table, sha256, guard


def verify_hash(path, expected):
    if not re.fullmatch(r"[0-9a-f]{64}", expected) or sha256(path) != expected:
        raise ValueError("Artifact/checkpoint SHA256 mismatch")


def validate_artifact(directory, job):
    directory = Path(directory).resolve()
    manifest = json.loads((directory / "run_manifest.json").read_text())
    if job["licence_status"] != "approved":
        raise ValueError("Method excluded: unresolved source licence")
    for k in ["method", "target", "strategy", "source_repository", "source_commit", "configuration_hash", "input_hashes", "context"]:
        if manifest.get(k) != job[k]:
            raise ValueError(f"Remote provenance mismatch: {k}")
    if manifest.get("assay_values_used") is not False or manifest.get("status") != "complete":
        raise ValueError("Assay leakage declaration or run status invalid")
    required = ["package_versions", "cuda_version", "gpu_model", "gpu_memory_bytes", "peak_allocated_gpu_memory_bytes", "start_time", "end_time", "elapsed_seconds", "checkpoint_hash", "training_seed", "generation_seeds", "generation_statistics", "output_hashes"]
    if any(k not in manifest for k in required):
        raise ValueError("Incomplete remote run manifest")
    if manifest["generation_seeds"] != job["generation_seeds"] or manifest["training_seed"] != job["training_seed"]:
        raise ValueError("Seed mismatch")
    if not re.fullmatch(r"[0-9a-f]{64}", manifest["checkpoint_hash"]):
        raise ValueError("Invalid checkpoint content hash")
    expected = {f"generations_{s}.tsv.gz" for s in job["generation_seeds"]} | {"heldout_scores.tsv"}
    if set(manifest["output_hashes"]) != expected:
        raise ValueError("Unexpected or missing returned files")
    for name, digest in manifest["output_hashes"].items():
        original = directory / name
        path = original.resolve()
        if not path.is_relative_to(directory) or original.is_symlink():
            raise ValueError("Artifact path escape")
        verify_hash(path, digest)
    train = set(job.get("training_sequences", []))
    for seed in job["generation_seeds"]:
        sequences = [r["sequence"] for r in rows(directory / f"generations_{seed}.tsv.gz")]
        if len(sequences) != job["generation_budget"] or len(set(sequences)) != len(sequences):
            raise ValueError("Unique generation budget not met")
        if any(set(s)-set("ACGT") or not job["minimum_length"] <= len(s) <= job["maximum_length"] or any(p in s for p in job["forbidden_primers"]) for s in sequences):
            raise ValueError("Invalid generated sequence or primer contamination")
    return manifest


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--directory", type=Path, required=True)
    p.add_argument("--job", type=Path, required=True)
    a = p.parse_args()
    job = json.loads(a.job.read_text())
    registered = list(rows(ROOT / "results/gpu_jobs.tsv"))
    if not any(r["method"] == job["method"] and r["target"] == job["target"] and r["strategy"] == job["strategy"] and r["job_sha256"] == sha256(a.job) for r in registered):
        raise ValueError("Job file differs from registered local export")
    # Validate the local bundle too; a copied training file must not drift after export.
    for name, digest in job["input_hashes"].items():
        verify_hash(a.job.parent / name, digest)
    manifest = validate_artifact(a.directory, job)
    expected_test = {r["sequence"] for r in rows(a.job.parent / "test_sequences.tsv")}
    scores = list(rows(a.directory / "heldout_scores.tsv"))
    if {r["sequence"] for r in scores} != expected_test or len(scores) != len(expected_test):
        raise ValueError("Held-out scores omit or duplicate test sequences")
    import math
    if any(not math.isfinite(float(r["score"])) or r["score_type"] != "posterior_mean_conditional_PHMM_log_probability_per_nt" for r in scores):
        raise ValueError("Unregistered or nonfinite model score")
    guard(sum((a.directory/n).stat().st_size for n in manifest["output_hashes"]))
    dest = ROOT / f'data/remote/{job["method"]}_{job["target"]}_{job["strategy"]}'
    dest.mkdir(parents=True, exist_ok=True)
    for name in list(manifest["output_hashes"]) + ["run_manifest.json"]:
        shutil.copy2(a.directory/name, dest/name)
    table(ROOT / f'results/import_{job["method"]}_{job["target"]}_{job["strategy"]}.tsv', [dict(method=job["method"],target=job["target"],strategy=job["strategy"],manifest_sha256=sha256(dest/"run_manifest.json"), status="validated", checkpoint_hash=manifest["checkpoint_hash"], note="Checkpoint hash is remote-recorded; weights intentionally not returned")])
    print("Imported validated model outputs; generated sequences remain candidates")


if __name__ == "__main__":
    main()
