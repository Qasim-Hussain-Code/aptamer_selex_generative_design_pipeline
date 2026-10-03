"""Report primary-source identity, repository licences and accession contradictions."""
import argparse
import csv
import json
import statistics
import xml.etree.ElementTree as ET
from collections import Counter
from common import ROOT, rows, table, now, sha256, failure, completed, stamp
from sequences import normalize
import importlib


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--offline", action="store_true", help="Validate cached evidence without network")
    p.parse_args()
    inputs = list((ROOT/"data/external").glob("*.xml")) + list((ROOT/"data/external").glob("*/repository.json")) + list((ROOT/"data/external").glob("*/tree.json")) + list((ROOT/"data/external").glob("*/commits.json")) + list((ROOT/"data/external/aptadiff/data/raw_data").glob("*")) + [ROOT/"results/local_software.tsv",ROOT/"config/datasets.tsv"]
    outputs=[ROOT/f"results/{n}.tsv" for n in ["citations","upstream_software","optional_method_exclusions","software_manifest","data_provenance","unresolved_items","data_terms"]]
    if completed("sources",inputs,outputs):return
    refs = [dict(method="RaptRanker", doi="10.1093/nar/gkaa484", title="RaptRanker: in silico RNA aptamer selection from HT-SELEX experiment based on local sequence and structure information", source="https://academic.oup.com/nar/article/48/14/e82/5855635", verification="primary fullTextXML and supplementary PDF", checked_date="2026-10-03"),
            dict(method="RaptGen", doi="10.1038/s43588-022-00249-6", title="Generative aptamer discovery using RaptGen", source="https://www.nature.com/articles/s43588-022-00249-6", verification="primary publisher article inspected", checked_date="2026-10-03"),
            dict(method="AptaDiff", doi="10.1093/bib/bbae517", title="AptaDiff: de novo design and optimization of aptamers based on diffusion models", source="https://pmc.ncbi.nlm.nih.gov/articles/PMC11491854/", verification="primary fullTextXML and publisher article", checked_date="2026-10-03"),
            dict(method="InstructNA", doi="10.1038/s43588-026-00965-3", title="De novo design of functional nucleic acids of aptamers", source="https://www.nature.com/articles/s43588-026-00965-3", verification="primary publisher article inspected; optional code excluded", checked_date="2026-10-03")]
    for name, doi in [("raptranker", refs[0]["doi"]), ("aptadiff", refs[2]["doi"])]:
        doc = ET.parse(ROOT / f"data/external/{name}.xml")
        ids = [x.text for x in doc.iter("article-id") if x.attrib.get("pub-id-type") == "doi"]
        if doi not in ids:
            raise ValueError(f"Primary XML DOI mismatch for {name}")
    table(ROOT / "results/citations.tsv", refs)
    software, exclusions = [], []
    for name in ["raptgen", "aptadiff", "instructna"]:
        d = ROOT / f"data/external/{name}"
        repo = json.loads((d / "repository.json").read_text())
        commit = json.loads((d / "commits.json").read_text())[0]["sha"]
        tree = json.loads((d / "tree.json").read_text())["tree"]
        paths = [x["path"] for x in tree if x["type"] == "blob"]
        licences = [x for x in paths if "licen" in x.lower() or "copying" in x.lower()]
        licence = repo["license"]["spdx_id"] if repo.get("license") and licences else "unclear"
        notes = "Licence file inspected; separately installed upstream code" if licence != "unclear" else "No explicit repository software licence at inspected commit"
        if name == "aptadiff":
            notes += "; paper Code availability declares MIT but repository contains no licence; excluded pending clarification"
        software.append(dict(tool=name, role="generator" if name != "instructna" else "optional generator", version="source commit",
             source=repo["html_url"], source_commit=commit, licence=licence, licence_url_or_source=repo["html_url"] + "/tree/" + commit,
             licence_checked_date="2026-10-03", execution="remote", redistributable="yes" if licence != "unclear" else "unclear", notes=notes))
        if licence == "unclear":
            exclusions.append(dict(method=name, reason=notes, source_commit=commit, checked_date="2026-10-03", status="excluded_no_code_or_weights_vendored"))
    table(ROOT / "results/upstream_software.tsv", software)
    table(ROOT / "results/optional_method_exclusions.tsv", exclusions)
    local = list(rows(ROOT / "results/local_software.tsv")) if (ROOT / "results/local_software.tsv").exists() else []
    table(ROOT / "results/software_manifest.tsv", local + software)
    provenance = []
    for target, acc in [("tg2", "DRA009383"), ("integrin", "DRA009384")]:
        found = list(rows(ROOT / f"data/external/{acc}.ena.tsv"))
        if not found:
            raise ValueError("Empty provider metadata")
        provenance.append(dict(file=acc, target=target, origin="DDBJ-submitted metadata via INSDC ENA mirror",
             evidence=found[0]["study_title"], accession_mapping="verified_primary_paper_and_provider_study_title", inconsistency="none",
             included="yes", reason="official primary benchmark input", minimum_length="", maximum_length="", primer_match_fraction=""))
    trim = importlib.import_module("05_preprocess_selex").trim
    fastq = importlib.import_module("05_preprocess_selex").fastq
    ds = {d["target"]: d for d in rows(ROOT / "config/datasets.tsv")}
    for letter, target in [("C", "tg2"), ("D", "integrin")]:
        path = ROOT / f"data/external/aptadiff/data/raw_data/dataset{letter}.fastq"
        lengths, passed, total = [], 0, 0
        for s in fastq(path):
            total += 1
            seq, reason = trim(s, ds[target])
            if seq:
                passed += 1
                lengths.append(len(seq))
        provenance.append(dict(file=f"AptaDiff/data/raw_data/dataset{letter}.fastq", target=target, origin="AptaDiff upstream repository",
             evidence="paper Table 1 maps C=TGM2 and D=integrin; measured primers agree", accession_mapping="consistent_library; round_and_subset_selection_unresolved",
             inconsistency="availability text assigns DRA accessions to A/B although Table 1 uses C/D", included="no",
             reason="use official full ENA rounds instead of undocumented repository subsets", minimum_length=min(lengths), maximum_length=max(lengths), primer_match_fraction=passed/total))
    for letter, target in [("A", "IGFBP3"), ("B", "PTK7")]:
        path = ROOT / f"data/external/aptadiff/data/raw_data/dataset{letter}_{target}_P6.csv"
        with open(path) as f:
            lens = [len(r["seq"]) for r in csv.DictReader(f)]
        provenance.append(dict(file=f"AptaDiff/data/raw_data/{path.name}", target=target, origin="authors' experimental library",
             evidence="paper SPR methods/Table 1 and repository filename; distinct primer designs", accession_mapping="not_DRA009383_or_DRA009384",
             inconsistency="unresolved public accession for A/B; paper nominal 36 nt compared with observed lengths",
             included="no", reason="distinct biological targets; no accession identity assumed", minimum_length=min(lens), maximum_length=max(lens), primer_match_fraction="variable_region_only"))
    table(ROOT / "results/data_provenance.tsv", provenance)
    table(ROOT / "results/unresolved_items.tsv", [dict(item="AptaDiff accession sentence", reason="A/B refer to IGFBP3/PTK7; DRA records identify TG2/integrin; paper Table 1 labels the public data C/D", action="Keep targets distinct; exclude ambiguous accession mapping"),
         dict(item="AptaDiff licence", reason="Paper states MIT; repository has no licence file", action="No upstream source used, modified or redistributed; remote job refuses"),
         dict(item="InstructNA licence and weights", reason="No repository licence at inspected commit; weight terms and footprint not established", action="Excluded before code use"),
         dict(item="DDBJ direct metadata endpoint", reason="Resource .xml route returned HTML; guessed FTP route returned 404", action="Use DDBJ-submitted study/experiment records from INSDC ENA mirror; record failures"),
         dict(item="RaptGen remote execution", reason="No hosted GPU run or CUDA compatibility test performed", action="Provide pinned job and validate returned outputs; no neural result claimed")])
    table(ROOT / "results/data_terms.tsv", [dict(component="RaptRanker article and supplement", terms="CC-BY-NC-4.0", source="https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7641312/fullTextXML", redistribution="primary PDF not redistributed; factual derived assay fields attributed; terms retained", checked_date="2026-10-03"),
         dict(component="DDBJ/ENA reads", terms="INSDC public archival data; no repository-specific licence inferred", source="https://www.ebi.ac.uk/ena/browser/about/terms-of-use", redistribution="raw and compact reads ignored; regenerate from accession manifests", checked_date="2026-10-03"),
         dict(component="AptaDiff repository datasets", terms="unclear separate dataset terms", source="https://github.com/wz-create/AptaDiff", redistribution="not redistributed; inspection summaries only", checked_date="2026-10-03")])
    print("Primary references, licence exclusions and biological provenance recorded")
    stamp("sources",inputs,outputs)


if __name__ == "__main__":
    main()
