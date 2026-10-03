import argparse
from common import ROOT, table


def main():
    p = argparse.ArgumentParser(description="Record why the optional tertiary arm has not run")
    p.parse_args()
    table(ROOT / "results/tertiary_preflight.tsv", [dict(method="none_selected", code_licence="not_evaluated", weight_terms="not_evaluated", hardware="no_hosted_run",
          held_out_endpoint="paired_delta_average_precision_on_fixed_assay_test_subset", proceed="no", reason="Optional arm disabled; no satisfactory method/licence/hardware preflight; no tertiary score produced")])
    print("Optional tertiary arm disabled by preregistered configuration")


if __name__ == "__main__":
    main()
