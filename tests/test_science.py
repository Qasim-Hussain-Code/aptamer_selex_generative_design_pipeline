import importlib
import math
from pathlib import Path
import pytest
import numpy as np
from common import rows
from sequences import families, assert_no_leakage, identity, deterministic_subset, reverse_complement
from model_utils import enrichment, assay_strength, binary_label

pre = importlib.import_module("05_preprocess_selex")


@pytest.fixture
def design():
    return dict(target="synthetic",forward_primer="ACGTTGCA",reverse_primer="GATCCAGT",fold_forward="ACGTTGCA",fold_reverse="GATCCAGT",minimum_length="7",maximum_length="9")


@pytest.mark.parametrize("core,reason",[("AAAACCCC","forward"),("AAAAAAA","forward"),("AAAAAAAAA","forward"),("AAAAAA","length_outside_published_range"),("AAAAAAAAAA","length_outside_published_range"),("AAAANCCC","ambiguous_base")])
def test_primer_length_and_ambiguity(design,core,reason):
    s,status=pre.trim(design["forward_primer"]+core+design["reverse_primer"],design)
    assert status==reason
    assert s==core if reason=="forward" else s is None


def test_wrong_primer(design):
    assert pre.trim("AAAAAAAAAAAACCCC"+design["reverse_primer"],design)==(None,"primer_mismatch")


def test_reverse_orientation(design):
    read=design["forward_primer"]+"AAAACCCC"+design["reverse_primer"]
    assert pre.trim(reverse_complement(read),design)==("AAAACCCC","reverse")


def test_counting_accountability(tmp_path,design):
    read=design["forward_primer"]+"AAAACCCC"+design["reverse_primer"]
    path=tmp_path/"tiny.fastq"
    path.write_text("".join(f"@{i}\n{s}\n+\n{'I'*len(s)}\n" for i,s in enumerate([read,read,"AAAA"])))
    result=pre.preprocess(path,tmp_path/"counts.tsv.gz",design,0,tmp_path/"qc.tsv")
    assert result["retained_reads"]==2 and result["discarded_reads"]==1
    assert result["raw_reads"]==result["retained_reads"]+result["discarded_reads"]
    assert next(rows(tmp_path/"counts.tsv.gz"))["read_count"]=="2"


def test_truncated_fastq_rejected(tmp_path):
    path=tmp_path/"bad.fastq"
    path.write_text("@bad\nACGT\n+\nII\n")
    with pytest.raises(ValueError):
        list(pre.fastq(path))


def test_enrichment_formula():
    assert enrichment(10,100,5,100,10,0.5)==pytest.approx(math.log2(10.5/5.5))


def test_zero_counts_are_finite():
    assert math.isfinite(enrichment(0,100,0,200,10,0.5))
    assert enrichment(5,100,0,100,10,0.5)>0


@pytest.mark.parametrize("late,early,vocabulary",[(0,1,10),(1,0,10),(1,1,0)])
def test_invalid_enrichment_depth(late,early,vocabulary):
    with pytest.raises(ValueError):
        enrichment(0,late,0,early,vocabulary)


def test_indels_and_transitive_families():
    seq=["AAAACCCC","AAAACCCCA","TAAACCCCA","GGGGTTTT"]
    assignment=families(seq,0.8)
    assert assignment[seq[0]]==assignment[seq[1]]==assignment[seq[2]]
    assert assignment[seq[0]]!=assignment[seq[3]]
    assert identity(seq[0],seq[1])==pytest.approx(8/9)
    assert families(reversed(seq),0.8)==assignment


def test_no_heldout_family_in_training():
    a=families(["AAAACCCC","AAAACCCA","GGGGTTTT"],0.8)
    assert_no_leakage(["GGGGTTTT"],["AAAACCCC"],a)
    with pytest.raises(AssertionError):
        assert_no_leakage(["AAAACCCA"],["AAAACCCC"],a)
    with pytest.raises(AssertionError):
        assert_no_leakage(["AAAACCCC"],["AAAACCCC"],a)


def test_kd_direction():
    assert assay_strength(1,"KD")>assay_strength(10,"KD")
    with pytest.raises(ValueError):
        assay_strength(0,"KD")


def test_response_direction():
    assert assay_strength(100,"SPR_end_of_injection_response")>assay_strength(-5,"SPR_end_of_injection_response")
    with pytest.raises(ValueError):
        assay_strength(1,"unknown")


def test_published_class_direction():
    assert binary_label("1")==1 and binary_label("0")==0
    with pytest.raises(ValueError):
        binary_label("binder")


def test_hash_subsampling_is_order_independent():
    seq=["AAAA","CCCC","GGGG","TTTT","ACGT"]
    assert deterministic_subset(seq,3,1729)==deterministic_subset(reversed(seq),3,1729)
    assert len(set(deterministic_subset(seq,3,1729)))==3


def test_seeded_generator_is_deterministic():
    from generation import fit_markov,generate
    model=fit_markov(["AAAACCCC","GGGGTTTT"],{"AAAACCCC":5,"GGGGTTTT":2})
    first=generate(model,11,30,10000)
    assert first==generate(model,11,30,10000)
    assert first!=generate(model,23,30,10000)


def test_tied_scores_use_prevalence():
    evaluate=importlib.import_module("15_evaluate")
    y=np.array([0,1,0,1])
    result=evaluate.metrics(y,np.zeros(4),np.array([0,10,2,4]),2)
    assert result["average_precision"]==0.5
    assert result["precision_top_k"]==0.5


def test_synthetic_end_to_end(tmp_path):
    from smoke import smoke
    assert smoke(tmp_path)["status"]=="passed"


def test_score_schema():
    scorer=importlib.import_module("12_score_candidates")
    row=dict(zip(scorer.SCORE_FIELDS,["tg2","0.8","family","markov","ACGT","fam1","-1","per_nt_log_probability","1","30","RU","harmonized_benchmark","0","10"]))
    scorer.validate_scores([row])
    with pytest.raises(ValueError):
        scorer.validate_scores([dict(row,score="nan")])
    with pytest.raises(ValueError):
        scorer.validate_scores([row,row])


@pytest.mark.parametrize("seed",range(10))
def test_bootstrap_ap_matches_sklearn(seed):
    from sklearn.metrics import average_precision_score
    rng=np.random.default_rng(seed)
    y=np.r_[0,1,rng.integers(0,2,38)]
    scores=rng.integers(0,5,40).astype(float)
    assert importlib.import_module("15_evaluate").average_precision(y,scores)==pytest.approx(average_precision_score(y,scores))
