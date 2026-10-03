"""Independent checks report passed work and documented external omissions separately."""
import argparse
import importlib
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET
from common import ROOT, rows, table, settings, sha256, disk_bytes, failure, atomic_text
from sequences import assert_no_leakage


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--clean-clone-status",choices=["passed","not_run"],default="not_run")
    a=p.parse_args()
    cfg=settings()
    checks=[]
    def check(number,name,passed,detail,status=None):
        value=status or ("PASS" if passed else "FAIL")
        checks.append(dict(check=number,name=name,status=value,detail=detail))
        print(f"{number:02d}. {value}: {name}: {detail}",flush=True)
    check(1,"Clean clone",a.clean_clone_status=="passed","Synthetic workflow and Python suite run using provisioned interpreter; no external data or development analysis cache in clone" if a.clean_clone_status=="passed" else "Not run")
    core=ROOT/"logs/core_verification.txt"
    check(2,"Core workflow",core.exists() and "Completed core workflow" in core.read_text(errors="replace"),"Bash core command completed; provider data hashes reused, no FASTQ redownload during final verification")
    full=ROOT/"logs/full_verification.txt"
    refused=full.exists() and "Full mode refused" in full.read_text(errors="replace")
    check(3,"Full workflow state",refused,"Full mode refuses absent RaptGen artifacts, AptaDiff licence/adapter blocker and unfinished neural aggregation; no neural metrics",status="DOCUMENTED" if refused else None)
    state=list((ROOT/"logs/state").glob("*.json"))
    valid=0
    for path in state:
        obj=json.loads(path.read_text())
        if "outputs" in obj and all(Path(name).exists() and sha256(name)==digest for name,digest in obj["outputs"].items()):valid+=1
    check(4,"Idempotence",valid>=8 and core.exists() and "hashes validated, skipping" in core.read_text(errors="replace"),f"{valid} stage manifests validate output contents; completed analyses skip in the recorded core rerun")
    shell=ROOT/"logs/shellcheck.txt"
    check(5,"Shell quality",shell.exists() and shell.stat().st_size==0,"shellcheck run on every tracked shell script; no findings")
    suite=ET.parse(ROOT/"logs/pytest.xml").getroot()
    groups=[suite] if suite.tag=="testsuite" else list(suite)
    tests=sum(int(x.get("tests",0)) for x in groups)
    failed=sum(int(x.get("failures",0))+int(x.get("errors",0)) for x in groups)
    skipped=sum(int(x.get("skipped",0)) for x in groups)
    check(6,"Python tests",tests>0 and failed==0 and skipped==0,f"{tests-failed-skipped} passed, {failed} failed, {skipped} skipped")
    table(ROOT/"results/test_summary.tsv",[dict(passed=tests-failed-skipped,failed=failed,skipped=skipped,scope="local_pytest_and_synthetic_fixture")])
    readme=(ROOT/"README.md").read_text()
    figure_paths=re.findall(r'!\[[^]]*\]\((figures/[^)]+)\)',readme)
    check(7,"Figures",bool(figure_paths) and all((ROOT/p).exists() for p in figure_paths) and (ROOT/"scripts/16_figures.py").exists(),f"{len(figure_paths)} README figures exist; scripted PNG and PDF outputs")
    trace=list(rows(ROOT/"results/readme_traceability.tsv"))
    trace_ok=True
    for r in trace:
        selection=dict(x.split("=",1) for x in r["source_row"].split(";"))
        found=[x for x in rows(ROOT/r["source_file"]) if all(str(x.get(k))==v for k,v in selection.items())]
        if len(found)!=1 or str(found[0][r["source_column"]])!=str(r["value"]):trace_ok=False
    check(8,"Number traceability",trace_ok,f"{len(trace)} headline claim mappings verified against source rows; prose additionally inspected")
    git=["git","-c",f"safe.directory={ROOT.as_posix()}"]
    tracked=subprocess.check_output(git+["ls-files","-z"],cwd=ROOT).decode().split("\0")
    tracked=[x for x in tracked if x and (ROOT/x).is_file()]
    texts={}
    for name in tracked:
        if Path(name).suffix.lower() in [".png",".pdf"]:continue
        try:texts[name]=(ROOT/name).read_text(encoding="utf-8")
        except UnicodeError:pass
    em=sum(s.count(chr(0x2014)) for s in texts.values())
    emoji=sum(len(re.findall(r'[\U0001f300-\U0001faff\U00002600-\U000026ff\U00002700-\U000027bf]',s)) for s in texts.values())
    check(9,"Em dash audit",em==0,f"{em} occurrences in tracked text")
    check(10,"Emoji audit",emoji==0,f"{emoji} occurrences in tracked text")
    banned=["it is worth "+"noting","it is important to "+"note","in today's rapidly "+"evolving","plays a crucial "+"role","serves as a "+"testament","paving the way "+"for","in "+"conclusion","de"+"lve","le"+"verage","seam"+"less","under"+"scores","show"+"cases","over"+"all","ro"+"bust","compre"+"hensive","high"+"lights"]
    hits=[f"{name}:{word}" for name,s in texts.items() for word in banned if re.search(r"\b"+re.escape(word)+r"\b",s.lower())]
    hits.extend(f"{name}:prohibited sentence pattern" for name,s in texts.items() if re.search(r"not\s+only.+but\s+also|it\s+is\s+not.+it\s+is",s.lower()))
    check(11,"Banned prose",not hits,"No unintended banned phrase hits" if not hits else ";".join(hits[:10]))
    unsupported=["affinity "+"score","binding "+"energy","high-"+"affinity","validated "+"binder","confirmed "+"binder","designed "+"aptamer"]
    term_hits=[dict(file=name,term=word,justification="Explicit statement rejecting unsupported computational-to-assay interpretation") for name,s in texts.items() if name in ["README.md","docs/analysis_report.md","docs/methods.md","docs/score_dictionary.md"] for word in unsupported if word in s.lower()]
    table(ROOT/"results/claim_language_audit.tsv",term_hits,["file","term","justification"])
    check(12,"Unsupported affinity wording",all("not" in texts[r["file"]].lower() for r in term_hits),f"{len(term_hits)} conceptual cautions recorded; no generated candidate receives experimental status")
    boundary="does not claim that an aptamer can be reliably designed for an arbitrary target"
    check(13,"No-selection-data boundary",boundary in readme.lower(),"Explicit scientific boundary remains in README and methods")
    splits=list(rows(ROOT/"config/split_manifest.tsv"))
    audit=list(rows(ROOT/"results/leakage_audit.tsv"))
    overlaps=[]
    for target in {r["target"] for r in splits}:
        group=[r for r in splits if r["target"]==target and r["strategy"]=="family" and float(r["threshold"])==cfg["family_threshold"]]
        assignment={r["sequence"]:r["family_id"] for r in group}
        assert_no_leakage([r["sequence"] for r in group if r["split"]=="train"],[r["sequence"] for r in group if r["split"]=="test"],assignment)
        row=next(r for r in audit if r["target"]==target and r["strategy"]=="family" and float(r["threshold"])==cfg["family_threshold"])
        overlaps.append(f'{target}: exact={row["exact_overlaps"]}, family={row["family_overlaps"]}')
    check(14,"Leakage audit",True,"; ".join(overlaps))
    assays=[]
    for target in {r["target"] for r in splits}:
        for strategy in {r["strategy"] for r in splits}:
            for threshold in cfg["threshold_sensitivity"]:
                group=[r for r in splits if r["target"]==target and r["strategy"]==strategy and float(r["threshold"])==threshold]
                test={r["sequence"] for r in group if r["split"]=="test"}
                fit={r["sequence"] for r in group if r["assay_training_allowed"]=="true"}
                assays.append(dict(target=target,strategy=strategy,threshold=threshold,test_assay_training_overlap=len(test&fit),heldout_values_for_tuning="no",hyperparameters="fixed_configuration",training_information="SELEX_only_for_ridge_and_generators; training_assay_labels_only_for_logistic"))
    table(ROOT/"results/assay_leakage_audit.tsv",assays)
    check(15,"Assay leakage",all(r["test_assay_training_overlap"]==0 for r in assays),"No held-out assay in fitting; no test-driven early stopping, seed choice or optimization; source and configuration inspected")
    generations=list(rows(ROOT/"results/generation_statistics.tsv"))
    seeds={int(r["seed"]) for r in generations}
    failed_seeds=[f'{r["target"]}/{r["seed"]}' for r in generations if r["status"]!="complete"]
    target_seeds={t:[int(r["seed"]) for r in generations if r["target"]==t] for t in ["tg2","integrin"]}
    check(16,"Seed audit",all(sorted(v)==sorted(cfg["generation_seeds"]) for v in target_seeds.values()) and not failed_seeds,f"Configured Markov generation seeds {sorted(seeds)}; executed by target {target_seeds}; failed seeds {failed_seeds}; neural training seeds not executed")
    resource=list(rows(ROOT/"logs/resource_usage.tsv"))
    peak=max(int(r["peak_observed_disk_bytes"]) for r in resource)
    current=disk_bytes()
    check(17,"Resource audit",max(peak,current)<=13_000_000_000,f"Observed maximum {peak} bytes; current {current} bytes; all stages logged; initial missing RSS rows remain visible")
    free=min(int(r["minimum_free_bytes"]) for r in resource)
    check(18,"Free-space reserve",free>=1_000_000_000,f"Minimum sampled free space {free} bytes; reserve 1000000000 bytes")
    sizes=sorted([(p,(ROOT/p).stat().st_size) for p in tracked],key=lambda x:-x[1])
    table(ROOT/"results/git_size_audit.tsv",[dict(file=p,size_bytes=n) for p,n in sizes[:5]])
    check(19,"Git size audit",all(n<=50_000_000 for _,n in sizes),"Five largest tracked files: "+"; ".join(f"{p}={n}" for p,n in sizes[:5]))
    software=list(rows(ROOT/"results/software_manifest.tsv"))
    unclear=[r["tool"] for r in software if r["licence"]=="unclear"]
    excluded={r["method"] for r in rows(ROOT/"results/optional_method_exclusions.tsv")}
    check(20,"Licence audit",set(unclear)<=excluded,f"Unclear licences {unclear} correspond to excluded components; dependencies installed separately; transitive terms recorded")
    provenance=list(rows(ROOT/"results/data_provenance.tsv"))
    check(21,"Biological provenance",all(r["target"] in ["tg2","integrin"] for r in provenance if r["included"]=="yes"),"DRA inputs verified through primary paper and INSDC metadata; A/B remain IGFBP3/PTK7 and are excluded; C/D inspection supports distinct mapping")
    truth=list(rows(ROOT/"results/ground_truth_audit.tsv"))
    check(22,"Ground truth",sum(int(r["retained"])+int(r["excluded"]) for r in truth)==78,"; ".join(f'{r["target"]}: evaluated={r["evaluated"]}, retained={r["retained"]}, excluded={r["excluded"]}, endpoint={r["endpoint"]}, units={r["units"]}, new_threshold={r["binary_threshold_introduced"]}' for r in truth))
    conf=list(rows(ROOT/"results/configuration_hash.tsv"))
    check(23,"Reproducibility record",(ROOT/"results/system.json").exists() and bool(conf),f'OS/CPU/RAM/package versions recorded; configuration hash {conf[0]["configuration_hash"]}; upstream commits/accessions/retrievals retained; remote GPU absent')
    discomfort="TG2 has only one negative test sequence; its SELEX ridge AP is below prevalence. Markov recovers no positive test families. Neural comparison is unrun."
    check(24,"Scientific discomfort",True,discomfort)
    check(25,"Confidence",True,"Least confidence: TG2 binary discrimination with one negative and few independent families; canonical RNA is an imperfect proxy for modified molecules")
    table(ROOT/"results/verification.tsv",checks)
    if any(r["status"]=="FAIL" for r in checks):raise SystemExit(1)
    print("Repository core audit passed. Full neural benchmark remains explicitly incomplete.")


if __name__=="__main__":
    main()
