"""First-order Markov null; count weighting is defined before held-out evaluation."""
from collections import Counter
import json
import math
import random
import numpy as np
from rapidfuzz import process
from rapidfuzz.distance import Levenshtein
from common import ROOT, rows, table, settings, dataset, sha256, atomic_text, completed, stamp
import importlib


def fit_markov(sequences, counts):
    starts = np.ones(4)
    transitions = np.ones((4, 4))
    bases = "ACGT"
    lengths = Counter()
    for s in sequences:
        n = counts[s]
        starts[bases.index(s[0])] += n
        lengths[len(s)] += n
        for a, b in zip(s[:-1], s[1:]):
            transitions[bases.index(a), bases.index(b)] += n
    return dict(start=(starts/starts.sum()).tolist(), transition=(transitions/transitions.sum(axis=1, keepdims=True)).tolist(), lengths=dict(lengths))


def log_probability(s, model):
    bases = "ACGT"
    p = math.log(model["start"][bases.index(s[0])])
    for a, b in zip(s[:-1], s[1:]):
        p += math.log(model["transition"][bases.index(a)][bases.index(b)])
    return p / len(s)


def generate(model, seed, budget, maximum):
    rng = random.Random(seed)
    seen, output = set(), []
    lengths = sorted(int(k) for k in model["lengths"])
    weights = [model["lengths"].get(k, model["lengths"].get(str(k))) for k in lengths]
    draws = 0
    while len(output) < budget and draws < maximum:
        draws += 1
        length = rng.choices(lengths, weights=weights)[0]
        s = rng.choices("ACGT", weights=model["start"])[0]
        for _ in range(length-1):
            s += rng.choices("ACGT", weights=model["transition"]["ACGT".index(s[-1])])[0]
        if s not in seen:
            seen.add(s)
            output.append(s)
    return output, draws


def main():
    cfg = settings()
    inputs = [ROOT / "config/benchmark.yml", ROOT / "config/split_manifest.tsv"] + list((ROOT / "data/processed").glob("*/*.counts.tsv.gz"))
    outputs = [ROOT / "results/generation_statistics.tsv", ROOT / "results/generation_budget_curves.tsv", ROOT / "results/markov_predictions.tsv", ROOT / "results/generator_metadata.tsv", ROOT / "results/novelty_distribution.tsv"]
    for target in ["tg2", "integrin"]:
        outputs.append(ROOT / f"data/processed/{target}/markov.json")
        for seed in cfg["generation_seeds"]:
            outputs.extend([ROOT / f"data/generated/{target}/markov_{seed}.tsv.gz", ROOT / f"data/generated/{target}/markov_{seed}.novelty.tsv.gz"])
    if completed("generation", inputs, outputs):
        return
    split = list(rows(inputs[1]))
    truth = list(rows(ROOT / "config/ground_truth_manifest.tsv"))
    metadata, stats, curves, predictions, novelty = [], [], [], [], []
    for target in sorted({r["target"] for r in split}):
        counts = importlib.import_module("07_build_splits").load_counts(target)
        late = counts[max(counts)]
        d = dataset(target)
        for threshold in cfg["threshold_sensitivity"]:
            for strategy in ["family", "random", "naive_sequence_random"]:
                group = [r for r in split if r["target"] == target and float(r["threshold"]) == threshold and r["strategy"] == strategy]
                train = sorted(r["sequence"] for r in group if r["split"] == "train")
                testrows = [r for r in group if r["split"] == "test"]
                model = fit_markov(train, late)
                for r in testrows:
                    tr = next(x for x in truth if x["target"] == target and x["sequence"] == r["sequence"])
                    predictions.append(dict(target=target, threshold=threshold, strategy=strategy, method="markov", sequence=r["sequence"], family_id=r["family_id"],
                         score=log_probability(r["sequence"], model), score_type="per_nt_log_probability", experimental_label=tr["original_experimental_label"], assay_measurement=tr["original_assay_measurement"],
                         assay_units="RU", context="harmonized_benchmark", assay_training_count=0, training_sequences=len(train)))
                if threshold != cfg["family_threshold"] or strategy != "family":
                    continue
                modelpath = ROOT / f"data/processed/{target}/markov.json"
                atomic_text(modelpath, json.dumps(model, sort_keys=True))
                test = [r["sequence"] for r in testrows]
                testfam = {r["sequence"]: r["family_id"] for r in testrows}
                training_families = {r["sequence"]: r["family_id"] for r in group}
                positive_fams = {testfam[r["sequence"]] for r in truth if r["target"] == target and r["sequence"] in testfam and r["original_experimental_label"] == "1"}
                for seed in cfg["generation_seeds"]:
                    sequences, draws = generate(model, seed, cfg["generation_budget"], cfg["generation_max_draws"])
                    dest = ROOT / f"data/generated/{target}/markov_{seed}.tsv.gz"
                    table(dest, [dict(sequence=s, score=log_probability(s, model)) for s in sequences])
                    identities, distances, recovered, recovered_families = [], [], set(), set()
                    nearest_records = []
                    for i, s in enumerate(sequences, 1):
                        nearest, sim, _ = process.extractOne(s, train, scorer=Levenshtein.normalized_similarity)
                        dist = Levenshtein.distance(s, nearest)
                        identities.append(sim)
                        distances.append(dist)
                        if s in testfam:
                            recovered.add(s)
                        for tr in test:
                            if Levenshtein.normalized_similarity(s, tr) + 1e-12 >= threshold:
                                recovered_families.add(testfam[tr])
                        nearest_records.append(dict(sequence=s, nearest_training_sequence=nearest, maximum_identity=sim, edit_distance=dist,
                             exact_training_membership=str(s in late and s in train).lower(), nearest_training_family=training_families[nearest], score=log_probability(s, model)))
                        if i in cfg["generation_curve_budgets"]:
                            curves.append(dict(target=target, method="markov", seed=seed, budget=i, exact_assay_recovery=len(recovered), family_recovery=len(recovered_families),
                                 positive_family_recovery=len(recovered_families & positive_fams), test_families=len(set(testfam.values())), positive_test_families=len(positive_fams)))
                    table(ROOT / f"data/generated/{target}/markov_{seed}.novelty.tsv.gz", nearest_records)
                    contamination = sum(d["fold_forward"] in s or d["fold_reverse"] in s for s in sequences)
                    stats.append(dict(target=target, method="markov", seed=seed, raw_generations=draws, valid_generations=draws, unique_generations=len(sequences), duplicates=draws-len(sequences),
                         requested_unique_budget=cfg["generation_budget"], status="complete" if len(sequences) == cfg["generation_budget"] else "budget_unattained",
                         exact_test_recovery=len(recovered), test_family_recovery=len(recovered_families), positive_test_family_recovery=len(recovered_families & positive_fams),
                         constant_region_contamination=contamination, mean_training_identity=float(np.mean(identities)), mean_nearest_edit_distance=float(np.mean(distances)),
                         mean_length=float(np.mean([len(s) for s in sequences])), gc_fraction=sum(s.count("G")+s.count("C") for s in sequences)/sum(map(len,sequences)),
                         exact_training_copies=sum(s in train for s in sequences), output_sha256=sha256(dest)))
                    for name, values in [("maximum_training_identity", identities), ("nearest_edit_distance", distances), ("length", list(map(len,sequences)))]:
                        for q in [0, 0.25, 0.5, 0.75, 1]:
                            novelty.append(dict(target=target, method="markov", seed=seed, quantity=name, quantile=q, value=float(np.quantile(values,q))))
                    metadata.append(dict(target=target, method="markov", seed=seed, training_seed=cfg["master_seed"], source="repository scripts/generation.py", source_commit="own_code", licence="MIT",
                         context="unsupervised_or_selex_only_generation", architecture="first_order_count_weighted_markov_add_one_smoothing", training_data_hash=sha256(inputs[1]),
                         checkpoint_hash=sha256(modelpath), environment="local core", hardware="CPU", training_sequences=len(train), test_assay_values_used="no", output_sha256=sha256(dest), status=stats[-1]["status"]))
                    print(f"Generated {target}/markov/{seed}: {len(sequences)} unique candidates", flush=True)
    for path, records in zip(outputs[:5], [stats, curves, predictions, metadata, novelty], strict=True):
        table(path, records)
    stamp("generation", inputs, outputs)
