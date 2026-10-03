import argparse
from common import ROOT, rows, table, completed, stamp
from generation import main as generate

SCORE_FIELDS = ["target", "threshold", "strategy", "method", "sequence", "family_id", "score", "score_type", "experimental_label", "assay_measurement", "assay_units", "context", "assay_training_count", "training_sequences"]


def validate_scores(records):
    import math
    seen = set()
    for r in records:
        if set(r) != set(SCORE_FIELDS) or not math.isfinite(float(r["score"])) or str(r["experimental_label"]) not in ["0", "1"]:
            raise ValueError("Score-table schema, direction or finite-value failure")
        key = tuple(r[k] for k in ["target", "threshold", "strategy", "method", "sequence"])
        if key in seen:
            raise ValueError("Duplicate score row")
        seen.add(key)


def main():
    p = argparse.ArgumentParser(description="Generate the null and standardize precisely named scores")
    p.parse_args()
    generate()
    from generation_qc import main as generation_qc
    generation_qc()
    inputs = [ROOT / f"results/{n}.tsv" for n in ["baseline_predictions", "structure_predictions", "markov_predictions"]]
    out = ROOT / "results/held_out_predictions.tsv"
    if completed("scores", inputs, [out]):
        return
    records = [r for path in inputs for r in rows(path)]
    validate_scores(records)
    table(out, records, SCORE_FIELDS)
    stamp("scores", inputs, [out])


if __name__ == "__main__":
    main()
