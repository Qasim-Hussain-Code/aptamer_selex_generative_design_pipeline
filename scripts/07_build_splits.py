"""Construct exact edit-distance components before inspecting assay values."""
import argparse
from collections import defaultdict
import hashlib
from pathlib import Path
from common import ROOT, rows, table, settings, completed, stamp, sha256, atomic_text
from sequences import families, deterministic_subset, assert_no_leakage


def load_counts(target, root=ROOT):
    rounds = {}
    for path in sorted((Path(root) / f"data/processed/{target}").glob("*.counts.tsv.gz")):
        rnd = int(path.name.split(".")[0])
        rounds[rnd] = {r["sequence"]: int(r["read_count"]) for r in rows(path)}
    if len(rounds) < 2:
        raise RuntimeError(f"Need at least two processed rounds for {target}")
    return rounds


def build(target, cfg, root=ROOT):
    root = Path(root)
    counts = load_counts(target, root)
    late = max(counts)
    eligible = [s for s, n in counts[late].items() if n >= 2]
    pool = deterministic_subset(eligible, int(cfg["subsampling"]["training_pool"]), cfg["subsampling"]["seed"])
    assays = [r for r in rows(root / "config/ground_truth_manifest.tsv") if r["target"] == target and r["inclusion_status"] == "included"]
    assay_sequences = {r["sequence"] for r in assays}
    assignment_primary = families(pool + list(assay_sequences), cfg["family_threshold"])
    assay_families = sorted({assignment_primary[s] for s in assay_sequences})
    # The number and identity of test families depend only on sequences, never on RU or labels.
    ntest = max(1, min(len(assay_families)-1, round(len(assay_families)*cfg["assay_test_family_fraction"])))
    test_families = set(sorted(assay_families, key=lambda f: hashlib.sha256(f'{cfg["master_seed"]}:{f}'.encode()).digest())[:ntest])
    test = {s for s in assay_sequences if assignment_primary[s] in test_families}
    records, audits = [], []
    for threshold in cfg["threshold_sensitivity"]:
        assign = assignment_primary if threshold == cfg["family_threshold"] else families(pool + list(assay_sequences), threshold)
        excluded_families = {assign[s] for s in test}
        for strategy in ["family", "random"]:
            train = {s for s in pool if s not in test and (strategy == "random" or assign[s] not in excluded_families)}
            if strategy == "family":
                assert_no_leakage(train, test, assign)
            for s in sorted(set(pool) | assay_sequences):
                split = "test" if s in test else "train" if s in train else "excluded"
                records.append(dict(target=target, threshold=threshold, strategy=strategy, sequence=s,
                     family_id=assign[s], split=split, experimental_label_status="measured" if s in assay_sequences else "unmeasured",
                     reason_for_assignment="fixed_primary_test_family" if split == "test" else "held_out_family_or_not_in_training_pool" if split == "excluded" else "eligible_late_round_training_sequence",
                     assay_training_allowed=str(s in assay_sequences and s not in test and (strategy == "random" or assign[s] not in excluded_families)).lower()))
            audits.append(dict(target=target, threshold=threshold, strategy=strategy,
                 training_sequences=len(train), test_sequences=len(test), test_families=len({assign[s] for s in test}),
                 exact_overlaps=len(train & test), family_overlaps=len({assign[s] for s in train} & {assign[s] for s in test}),
                 excluded_training_sequences=len(pool)-len(train), late_round_eligible=len(eligible), selected_pool=len(pool),
                 excluded_from_pool=len(eligible)-len(pool)))
    return records, audits


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", choices=["all", "tg2", "integrin"], default="all")
    a = p.parse_args()
    cfg = settings()
    targets = ["tg2", "integrin"] if a.target == "all" else [a.target]
    inputs = [ROOT / "config/benchmark.yml", ROOT / "config/ground_truth_manifest.tsv"] + list((ROOT / "data/processed").glob("*/*.counts.tsv.gz"))
    outputs = [ROOT / "config/split_manifest.tsv", ROOT / "results/leakage_audit.tsv"]
    if completed("splits", inputs, outputs):
        return
    recs, audits = [], []
    for t in targets:
        print(f"Constructing exact components for {t}", flush=True)
        r, au = build(t, cfg)
        recs.extend(r)
        audits.extend(au)
    table(outputs[0], recs)
    table(ROOT / "results/family_assignments.tsv", recs)
    table(outputs[1], audits)
    table(ROOT / "results/exclusions.tsv", [dict(target=r["target"], threshold=r["threshold"], strategy=r["strategy"], item=r["sequence"], reason=r["reason_for_assignment"]) for r in recs if r["split"] == "excluded"])
    stamp("splits", inputs, outputs)


if __name__ == "__main__":
    main()
