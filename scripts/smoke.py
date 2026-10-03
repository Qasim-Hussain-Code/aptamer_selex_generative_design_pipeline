"""Tiny synthetic end-to-end software check; no experimental claims."""
import argparse
import importlib
import json
from pathlib import Path
import tempfile
import numpy as np
from common import ROOT, rows, atomic_text
from sequences import families, assert_no_leakage
from generation import fit_markov, generate
from model_utils import enrichment


def smoke(directory):
    directory=Path(directory)
    d=dict(target="synthetic",forward_primer="ACGTTGCA",reverse_primer="GATCCAGT",minimum_length="7",maximum_length="9",fold_forward="ACGTTGCA",fold_reverse="GATCCAGT")
    seqs=["AAAACCCC","AAAACCCC","TTTTGGGG","TTTTGGGA","ACGTACGT"]
    reads=[d["forward_primer"]+s+d["reverse_primer"] for s in seqs]+["AAAA",d["forward_primer"]+"NNNNNNNN"+d["reverse_primer"]]
    fastq=directory/"fixture.fastq"
    fastq.write_text("".join(f"@synthetic{i}\n{s}\n+\n{'I'*len(s)}\n" for i,s in enumerate(reads)))
    qc=importlib.import_module("05_preprocess_selex").preprocess(fastq,directory/"counts.tsv.gz",d,1,directory/"qc.tsv")
    counted={r["sequence"]:int(r["read_count"]) for r in rows(directory/"counts.tsv.gz")}
    assert counted["AAAACCCC"]==2 and qc["raw_reads"]==7 and qc["retained_reads"]==5
    fam=families(list(counted),0.8)
    test={"TTTTGGGG"}
    train=[s for s in counted if fam[s] not in {fam[t] for t in test}]
    assert_no_leakage(train,test,fam)
    null=fit_markov(train,counted)
    candidates,draws=generate(null,11,20,10000)
    assert len(candidates)==len(set(candidates))==20
    structure=importlib.import_module("09_secondary_structure").fold("ACGTACGT","","",37)
    assert structure["base_pair_count"]>=0
    from sklearn.linear_model import Ridge
    from model_utils import kmer_features
    model=Ridge(alpha=10).fit(kmer_features(train),[enrichment(counted[s],5,0,5,4) for s in train])
    assert np.isfinite(model.predict(kmer_features(list(test)))).all()
    return dict(scope="synthetic_software_validation_only",raw_reads=7,retained_reads=5,unique_sequences=4,
                family_overlap=0,unique_candidates=20,folding_version=structure["viennarna_version"],status="passed")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    (ROOT/".cache").mkdir(exist_ok=True)
    with tempfile.TemporaryDirectory(dir=ROOT/".cache") as tmp:
        result=smoke(tmp)
    atomic_text(ROOT/"logs/smoke_summary.json",json.dumps(result,indent=2))
    print("Synthetic end-to-end smoke test passed; no public data or GPU required")


if __name__=="__main__":
    main()
