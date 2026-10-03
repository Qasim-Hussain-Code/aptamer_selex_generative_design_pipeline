"""ViennaRNA canonical-RNA proxy on the assay construct, excluding the poly(A) tether."""
import argparse
import RNA
from common import ROOT, rows, table, settings, dataset, completed, stamp


def fold(sequence, forward, reverse, temperature=37):
    # Constant regions can pair with the variable region; folding the core alone loses that context.
    construct = (forward + sequence + reverse).replace("T", "U")
    md = RNA.md()
    md.temperature = temperature
    fc = RNA.fold_compound(construct, md)
    structure, mfe = fc.mfe()
    fc.exp_params_rescale(mfe)
    fc.pf()
    count = structure.count("(")
    return dict(construct=construct, dot_bracket=structure, mfe_kcal_mol=mfe,
                base_pair_count=count, paired_fraction=2*count/len(construct),
                ensemble_diversity=fc.mean_bp_distance(), viennarna_version=RNA.__version__,
                model="canonical_RNA_secondary_structure_proxy_Turner2004", temperature_c=temperature)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    inputs = [ROOT / "config/split_manifest.tsv", ROOT / "config/datasets.tsv", ROOT / "config/benchmark.yml"]
    out = ROOT / "results/secondary_structure.tsv"
    if completed("secondary", inputs, [out]):
        return
    cfg = settings()
    selected = sorted({(r["target"], r["sequence"]) for r in rows(inputs[0])})
    results = []
    for i, (target, s) in enumerate(selected):
        d = dataset(target)
        results.append(dict(target=target, sequence=s, **fold(s, d["fold_forward"], d["fold_reverse"], cfg["folding_temperature_c"])))
        if i % 500 == 0:
            print(f"Folded {i}/{len(selected)} constructs", flush=True)
    table(out, results)
    stamp("secondary", inputs, [out])


if __name__ == "__main__":
    main()
