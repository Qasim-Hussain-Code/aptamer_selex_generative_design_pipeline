"""Extract Tables S4/S5 from the primary PDF, retaining published labels and RU."""
import argparse
import re
from pathlib import Path
from common import ROOT, rows, table, dataset, sha256, completed, stamp


def parse_text(text):
    pattern = r"(Data([12])-(\d+))\s+([ACGT]+)\s+(\d{6})\s+(-?\d+\.\d+)\s+([01])"
    found = re.findall(pattern, text)
    if len(found) != 78 or len({r[0] for r in found}) != 78:
        raise ValueError(f"Expected 78 distinct primary-source assay rows, parsed {len(found)}")
    return found


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--pdf", type=Path, default=ROOT / "data/external/raptranker_supplement.pdf")
    a = p.parse_args()
    out = ROOT / "config/ground_truth_manifest.tsv"
    inputs = [a.pdf, ROOT / "config/datasets.tsv"]
    audit = ROOT / "results/ground_truth_audit.tsv"
    outputs=[out,audit,ROOT / "results/experimental_ground_truth.tsv"]
    if completed("ground_truth", inputs, outputs):
        return
    from pypdf import PdfReader
    text = "\n".join(page.extract_text() for page in PdfReader(a.pdf).pages)
    records = []
    for ident, number, index, seq, date, ru, label in parse_text(text):
        target = "tg2" if number == "1" else "integrin"
        d = dataset(target)
        valid = int(d["minimum_length"]) <= len(seq) <= int(d["maximum_length"])
        records.append(dict(target=target, dataset_accession=d["accession"], assay_id=ident,
             sequence=seq, sequence_length=len(seq), original_assay_measurement=ru,
             assay_units="RU", original_experimental_label=label, endpoint="SPR_end_of_injection_response",
             direction="higher_is_stronger_response", source_table="S4" if number == "1" else "S5",
             source_citation="10.1093/nar/gkaa484", source_pdf_sha256=sha256(a.pdf),
             acquisition_date=date, notes="Published Positive/Negative (1/0); no new threshold. RU is not KD.",
             inclusion_status="included" if valid else "excluded", exclusion_reason="" if valid else "outside_published_length_range"))
    table(out, records)
    table(ROOT / "results/experimental_ground_truth.tsv", records)
    table(audit, [dict(target=t, evaluated=sum(r["target"] == t for r in records),
           retained=sum(r["target"] == t and r["inclusion_status"] == "included" for r in records),
           excluded=sum(r["target"] == t and r["inclusion_status"] != "included" for r in records),
           positives=sum(r["target"] == t and r["original_experimental_label"] == "1" for r in records),
           endpoint="published_binary_label_and_continuous_RU", units="RU", binary_threshold_introduced="no",
           original_label_criterion="response>25_RU_in_primary_paper_Table_1", exclusion_reasons="none") for t in ["tg2", "integrin"]])
    stamp("ground_truth", inputs, outputs)
    print(f"Reconstructed {len(records)} experimental rows directly from primary supplement")


if __name__ == "__main__":
    main()
