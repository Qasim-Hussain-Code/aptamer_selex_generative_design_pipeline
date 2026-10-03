"""Fit on training sequences only; observed assay-sequence counts are separate references."""
import argparse
import importlib
from pathlib import Path
import numpy as np
from sklearn.linear_model import Ridge, LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler
from common import ROOT, rows, table, settings, completed, stamp, failure
from model_utils import kmer_features, enrichment

load_counts = importlib.import_module("07_build_splits").load_counts


def predict_group(target, threshold, strategy, split, cfg, structures=None):
    rounds = load_counts(target)
    early, late = sorted(rounds)[-2:]
    ce, cl = rounds[early], rounds[late]
    ne, nl = sum(ce.values()), sum(cl.values())
    vocab = len(set(ce) | set(cl))
    train = sorted(r["sequence"] for r in split if r["split"] == "train")
    testrows = [r for r in split if r["split"] == "test"]
    test = [r["sequence"] for r in testrows]
    if not train or not test:
        raise RuntimeError(f"Empty training or test pool: {target}/{threshold}/{strategy}")
    truth = {r["sequence"]: r for r in rows(ROOT / "config/ground_truth_manifest.tsv") if r["target"] == target}
    x = kmer_features(train, cfg["kmer_sizes"])
    xt = kmer_features(test, cfg["kmer_sizes"])
    # Enrichment is a training response from SELEX. The ridge prediction has no assay units.
    y = np.asarray([enrichment(cl.get(s, 0), nl, ce.get(s, 0), ne, vocab, cfg["enrichment_pseudocount"]) for s in train])
    if structures is not None:
        features = cfg["structural_features"]
        x = np.column_stack([x, [[float(structures[(target,s)][f]) for f in features] for s in train]])
        xt = np.column_stack([xt, [[float(structures[(target,s)][f]) for f in features] for s in test]])
    model = make_pipeline(StandardScaler(), Ridge(alpha=cfg["ridge_alpha"]))
    model.fit(x, y)
    preds = {"kmer_structure_ridge" if structures is not None else "kmer_ridge": model.predict(xt)}
    if structures is None:
        # Exact test sequences have been removed from the harmonized training pool.
        # Keep the resulting constant baseline visible. Observed readouts use the full count table.
        preds["frequency"] = np.zeros(len(test))
        preds["enrichment"] = np.asarray([enrichment(0, nl, 0, ne, vocab, cfg["enrichment_pseudocount"])]*len(test))
        preds["frequency_observed"] = np.asarray([cl.get(s, 0)/nl for s in test])
        preds["enrichment_observed"] = np.asarray([enrichment(cl.get(s,0), nl, ce.get(s,0), ne, vocab, cfg["enrichment_pseudocount"]) for s in test])
    assay_train = sorted(r["sequence"] for r in split if r["assay_training_allowed"] == "true")
    if len(assay_train) >= 2 and len({truth[s]["original_experimental_label"] for s in assay_train}) == 2:
        xa = kmer_features(assay_train, cfg["kmer_sizes"])
        if structures is not None:
            xa = np.column_stack([xa, [[float(structures[(target,s)][f]) for f in cfg["structural_features"]] for s in assay_train]])
        classifier = make_pipeline(StandardScaler(), LogisticRegression(C=0.1, random_state=cfg["master_seed"], max_iter=2000))
        classifier.fit(xa, [int(truth[s]["original_experimental_label"]) for s in assay_train])
        preds["kmer_structure_logistic" if structures is not None else "kmer_logistic"] = classifier.decision_function(xt)
    else:
        failure("baselines", f"{target}/{threshold}/{strategy}/supervised", ValueError("Training assays lack both published classes"), "retain other methods and report exclusion", target, "excluded")
    records = []
    for method, scores in preds.items():
        for s, row, score in zip(test, testrows, scores, strict=True):
            tr = truth[s]
            records.append(dict(target=target, threshold=threshold, strategy=strategy, method=method,
                 sequence=s, family_id=row["family_id"], score=float(score),
                 score_type="decision_function" if "logistic" in method else "predicted_log2_enrichment" if "ridge" in method else "relative_frequency" if method.startswith("frequency") else "log2_enrichment",
                 experimental_label=tr["original_experimental_label"], assay_measurement=tr["original_assay_measurement"], assay_units=tr["assay_units"],
                 context="observational_reference" if method.endswith("observed") else "harmonized_benchmark",
                 assay_training_count=len(assay_train) if "logistic" in method else 0,
                 training_sequences=len(train)))
    return records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--target", choices=["all", "tg2", "integrin"], default="all")
    a = p.parse_args()
    inputs = [ROOT / "config/split_manifest.tsv", ROOT / "config/benchmark.yml", ROOT / "config/ground_truth_manifest.tsv"] + list((ROOT / "data/processed").glob("*/*.counts.tsv.gz"))
    out = ROOT / "results/baseline_predictions.tsv"
    if completed("baselines", inputs, [out]):
        return
    cfg = settings()
    splits = list(rows(inputs[0]))
    predictions = []
    for t in ["tg2", "integrin"] if a.target == "all" else [a.target]:
        for threshold in cfg["threshold_sensitivity"]:
            for strategy in ["family", "random"]:
                group = [r for r in splits if r["target"] == t and float(r["threshold"]) == threshold and r["strategy"] == strategy]
                predictions += predict_group(t, threshold, strategy, group, cfg)
    table(out, predictions)
    stamp("baselines", inputs, [out])


if __name__ == "__main__":
    main()
