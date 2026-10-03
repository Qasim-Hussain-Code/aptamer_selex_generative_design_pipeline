"""Family bootstrap intervals and paired differences use the same held-out sequences."""
import argparse
from collections import defaultdict
import warnings
import numpy as np
from scipy.stats import spearmanr, kendalltau
from sklearn.metrics import average_precision_score, roc_auc_score
from common import ROOT, rows, table, settings, completed, stamp


def metrics(y, score, response, top_k=10):
    ap = average_precision_score(y, score) if len(set(y)) > 1 else float("nan")
    auc = roc_auc_score(y, score) if len(set(y)) > 1 else float("nan")
    rho = float(spearmanr(response, score).statistic) if len(set(score)) > 1 and len(set(response)) > 1 else float("nan")
    tau = float(kendalltau(response, score).statistic) if len(set(score)) > 1 and len(set(response)) > 1 else float("nan")
    # Fractional inclusion at tied top-k boundaries avoids choosing an arbitrary sequence order.
    k = min(top_k, len(y))
    boundary = np.sort(score)[-k]
    above, tied = score > boundary, score == boundary
    precision = (y[above].sum() + (k-above.sum())*y[tied].mean())/k
    return dict(average_precision=float(ap), auroc=float(auc), spearman=rho, kendall=tau, precision_top_k=float(precision), prevalence=float(np.mean(y)))


def bootstrap_indices(families, replicates, seed):
    rng = np.random.default_rng(seed)
    unique = sorted(set(families))
    groups = {f: np.flatnonzero(np.asarray(families) == f) for f in unique}
    for _ in range(replicates):
        yield np.concatenate([groups[f] for f in rng.choice(unique, len(unique), replace=True)])


def ci(values):
    valid = np.asarray([v for v in values if np.isfinite(v)])
    return (float(np.quantile(valid, 0.025)), float(np.quantile(valid, 0.975)), len(valid)) if len(valid) else (float("nan"), float("nan"), 0)


def average_precision(y, score):
    """Noninterpolated AP with tied thresholds, checked against scikit-learn in tests."""
    if len(set(y))<2:
        return float("nan")
    order=np.argsort(-score,kind="stable")
    ranked=y[order]
    boundary=np.r_[np.flatnonzero(np.diff(score[order])),len(score)-1]
    positive=np.cumsum(ranked)[boundary]
    increments=np.diff(np.r_[0,positive])
    return float(np.sum(positive/(boundary+1)*increments)/positive[-1])


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    cfg = settings()
    inputs = [ROOT / "results/held_out_predictions.tsv", ROOT / "config/benchmark.yml", ROOT / "results/generation_statistics.tsv"]
    outputs = [ROOT / "results/evaluation_metrics.tsv", ROOT / "results/bootstrap_intervals.tsv", ROOT / "results/paired_differences.tsv", ROOT / "results/seed_variance.tsv"]
    if completed("evaluation", inputs, outputs):
        return
    grouped = defaultdict(list)
    for r in rows(inputs[0]):
        grouped[(r["target"], r["threshold"], r["strategy"], r["method"])].append(r)
    estimates, intervals, paired, arrays, boot = [], [], [], {}, {}
    warnings.filterwarnings("ignore", category=UserWarning)
    for key, records in grouped.items():
        records.sort(key=lambda r: r["sequence"])
        target, threshold, strategy, method = key
        y = np.asarray([int(r["experimental_label"]) for r in records])
        score = np.asarray([float(r["score"]) for r in records])
        response = np.asarray([float(r["assay_measurement"]) for r in records])
        fam = [r["family_id"] for r in records]
        arrays[key] = (y, score, response, fam, [r["sequence"] for r in records])
        values = []
        # All methods in a target/split use the same seeded family bootstrap samples.
        for idx in bootstrap_indices(fam, cfg["bootstrap_replicates"], cfg["master_seed"]):
            values.append(average_precision(y[idx], score[idx]))
        boot[key] = values
        point = metrics(y, score, response, cfg["top_k"])
        for metric, value in point.items():
            estimates.append(dict(target=target, threshold=threshold, strategy=strategy, method=method, metric=metric, value=value,
                 test_sequences=len(y), test_families=len(set(fam)), context=records[0]["context"], status="defined" if np.isfinite(value) else "undefined_constant_or_single_class"))
        lower, upper, valid = ci(values)
        intervals.append(dict(target=target, threshold=threshold, strategy=strategy, method=method, metric="average_precision", estimate=point["average_precision"],
             lower_95=lower, upper_95=upper, bootstrap_requested=cfg["bootstrap_replicates"], bootstrap_valid=valid,
             bootstrap_undefined=cfg["bootstrap_replicates"]-valid, resampling_unit="sequence_family", independent_families=len(set(fam))))
    for key in arrays:
        target, threshold, strategy, method = key
        comparisons = []
        if strategy == "random":
            comparisons.append(((target, threshold, "family", method), "random_minus_family"))
        if method in ["kmer_structure_ridge", "kmer_structure_logistic"]:
            comparisons.append(((target, threshold, strategy, method.replace("structure_", "")), "structure_minus_sequence"))
        for other, comparison in comparisons:
            if other not in arrays:
                continue
            left, right = arrays[key], arrays[other]
            if left[4] != right[4]:
                raise ValueError("Paired comparison uses different experimental sequences")
            # Use grouped-side families for both sides of the split contrast.
            fam = right[3] if comparison == "random_minus_family" else left[3]
            differences = []
            for idx in bootstrap_indices(fam, cfg["bootstrap_replicates"], cfg["master_seed"]):
                if len(set(left[0][idx])) < 2:
                    differences.append(float("nan"))
                else:
                    differences.append(average_precision(left[0][idx], left[1][idx])-average_precision(right[0][idx], right[1][idx]))
            lower, upper, valid = ci(differences)
            delta = average_precision_score(left[0], left[1])-average_precision_score(right[0], right[1]) if len(set(left[0])) > 1 else float("nan")
            paired.append(dict(target=target, threshold=threshold, strategy=strategy, method=method, comparison=comparison, metric="average_precision",
                 difference=delta, lower_95=lower, upper_95=upper, bootstrap_valid=valid, bootstrap_requested=cfg["bootstrap_replicates"], resampling_unit="sequence_family"))
    stats = list(rows(inputs[2]))
    variance = []
    for target in sorted({r["target"] for r in stats}):
        selected = [r for r in stats if r["target"] == target]
        for quantity in ["exact_test_recovery", "test_family_recovery", "positive_test_family_recovery", "mean_training_identity"]:
            vals = [float(r[quantity]) for r in selected]
            variance.append(dict(target=target, method="markov", quantity=quantity, configured_seeds=",".join(map(str,cfg["generation_seeds"])), executed_seeds=",".join(r["seed"] for r in selected),
                 failed_seeds=",".join(r["seed"] for r in selected if r["status"] != "complete"), mean=float(np.mean(vals)), standard_deviation=float(np.std(vals,ddof=1)),
                 minimum=min(vals), maximum=max(vals), variance_type="generation_seed_only; single_fixed_training_model"))
    for path, records in zip(outputs, [estimates, intervals, paired, variance], strict=True):
        table(path, records)
    stamp("evaluation", inputs, outputs)
    print(f"Evaluated {len(grouped)} method/split/threshold groups")


if __name__ == "__main__":
    main()
