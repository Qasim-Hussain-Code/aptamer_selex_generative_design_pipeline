"""Small I/O primitives shared by the numbered Bash-invoked stages."""
import csv
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import tempfile
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parents[1]


def now():
    return datetime.now(timezone.utc).isoformat()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def rows(path):
    path = Path(path)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8", newline="") as f:
        yield from csv.DictReader(f, delimiter="\t")


def table(path, data, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = list(data) if fields is None else data
    if fields is None:
        fields = list(data[0]) if data else []
    tmp = path.with_name(path.name + ".tmp")
    opener = gzip.open if path.suffix == ".gz" else open
    try:
        with opener(tmp, "wt", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fields, delimiter="\t", lineterminator="\n")
            w.writeheader()
            w.writerows(data)
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)


def atomic_text(path, text):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".tmp")
    tmp.write_text(text, encoding="utf-8", newline="\n")
    os.replace(tmp, path)


def append(path, record):
    path = Path(path)
    existing = list(rows(path)) if path.exists() else []
    table(path, existing + [record], list(record))


def failure(stage, item, error, action="stop", target="all", status="blocked"):
    append(ROOT / "logs/failures.tsv", dict(timestamp=now(), stage=stage,
           target=target, item=item, error_class=type(error).__name__,
           error=str(error)[:700], action_taken=action, final_status=status))


def download(url, dest, accession="", md5="", expected_bytes=None):
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_name(dest.name + ".part")
    req = urllib.request.Request(url, headers={"User-Agent": "aptamer-benchmark/1.0"})
    try:
        with urllib.request.urlopen(req, timeout=90) as response, open(tmp, "wb") as f:
            shutil.copyfileobj(response, f, 1024 * 1024)
        if expected_bytes is not None and tmp.stat().st_size != int(expected_bytes):
            raise ValueError("Downloaded byte count differs from provider manifest")
        local_md5 = hashlib.md5(tmp.read_bytes()).hexdigest() if tmp.stat().st_size < 10_000_000 else file_md5(tmp)
        if md5 and local_md5 != md5:
            raise ValueError("Provider MD5 mismatch")
        digest = sha256(tmp)
        os.replace(tmp, dest)
        append(ROOT / "results/download_manifest.tsv", dict(url=url, accession=accession,
               retrieval_timestamp=now(), size_bytes=dest.stat().st_size,
               provider_checksum=md5, local_md5=local_md5, local_sha256=digest,
               local_path=str(dest.relative_to(ROOT)) if dest.is_relative_to(ROOT) else str(dest)))
    finally:
        tmp.unlink(missing_ok=True)
    return dest


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def settings(root=ROOT):
    import yaml
    cfg = yaml.safe_load((Path(root) / "config/benchmark.yml").read_text())
    if "MASTER_SEED" in os.environ:
        cfg["master_seed"] = int(os.environ["MASTER_SEED"])
    return cfg


def dataset(target, root=ROOT):
    return next(r for r in rows(Path(root) / "config/datasets.tsv") if r["target"] == target)


def disk_bytes(root=ROOT):
    total = 0
    for parent, _, files in os.walk(root):
        for name in files:
            try:
                total += (Path(parent) / name).stat().st_size
            except FileNotFoundError:
                pass
    return total


def guard(extra=0, root=ROOT, safety=1.1):
    used = disk_bytes(root)
    free = shutil.disk_usage(root).free
    budget = int(os.environ.get("DISK_BUDGET_BYTES", 13_000_000_000))
    reserve = int(os.environ.get("RESERVE_BYTES", 1_000_000_000))
    projected = used + int(extra * safety)
    if projected > budget or free - int(extra * safety) < reserve:
        shortfall = max(projected - budget, reserve - (free - int(extra * safety)), 0)
        raise RuntimeError(f"Resource refusal: free={free} budget={budget} projected_peak={projected} "
                           f"safety_margin={safety - 1:.2f} shortfall={shortfall}. "
                           "Disable optional tertiary/checkpoint retention. A smaller preregistered "
                           "training pool may fit; full FASTQ download still requires one raw file.")
    return projected


STAGE_CODE = {
    "sources": ["01_verify_sources.py"],
    "ground_truth": ["06_build_ground_truth.py"],
    "splits": ["07_build_splits.py", "sequences.py"],
    "baselines": ["08_baselines.py", "07_build_splits.py", "model_utils.py"],
    "secondary": ["09_secondary_structure.py"],
    "ablation": ["13_structure_ablation.py", "08_baselines.py", "model_utils.py"],
    "generation": ["generation.py", "07_build_splits.py", "sequences.py"],
    "scores": ["12_score_candidates.py"],
    "evaluation": ["15_evaluate.py"],
}


def code_hashes(stage, root=ROOT):
    names = ["common.py"] + STAGE_CODE.get(stage, [])
    return {name: sha256(Path(root) / "scripts" / name) for name in names}


def stamp(stage, inputs, outputs, root=ROOT):
    record = dict(inputs={str(p): sha256(p) for p in inputs},
                  outputs={str(p): sha256(p) for p in outputs}, configuration=settings(root),
                  code=code_hashes(stage, root), time=now())
    atomic_text(Path(root) / f"logs/state/{stage}.json", json.dumps(record, indent=2))


def completed(stage, inputs, outputs, root=ROOT):
    path = Path(root) / f"logs/state/{stage}.json"
    if not path.exists():
        return False
    try:
        state = json.loads(path.read_text())
        valid = state["inputs"] == {str(p): sha256(p) for p in inputs} and state["outputs"] == {str(p): sha256(p) for p in outputs} and state.get("configuration") == settings(root) and state.get("code") == code_hashes(stage, root)
    except (OSError, ValueError, KeyError):
        return False
    if valid:
        print(f"{stage}: hashes validated, skipping")
    return valid
