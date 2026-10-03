import argparse
import importlib
from common import ROOT, rows, table, settings, completed, stamp


def main():
    p = argparse.ArgumentParser(description="Fit identical splits with additional canonical-RNA features")
    p.parse_args()
    inputs = [ROOT / "config/split_manifest.tsv", ROOT / "config/benchmark.yml", ROOT / "config/ground_truth_manifest.tsv", ROOT / "results/secondary_structure.tsv"]
    out = ROOT / "results/structure_predictions.tsv"
    if completed("ablation", inputs, [out]):
        return
    cfg = settings()
    splits = list(rows(inputs[0]))
    structures = {(r["target"], r["sequence"]): r for r in rows(inputs[3])}
    predictions = []
    for t in sorted({r["target"] for r in splits}):
        for threshold in cfg["threshold_sensitivity"]:
            for strategy in ["family", "random", "naive_sequence_random"]:
                group = [r for r in splits if r["target"] == t and float(r["threshold"]) == threshold and r["strategy"] == strategy]
                predictions += importlib.import_module("08_baselines").predict_group(t, threshold, strategy, group, cfg, structures)
    table(out, predictions)
    stamp("ablation", inputs, [out])


if __name__ == "__main__":
    main()
