"""Additional composition and diversity summaries, independent of experimental labels."""
import argparse
import random
from collections import Counter
from rapidfuzz.distance import Levenshtein
from common import ROOT, rows, table, settings


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.parse_args()
    cfg=settings()
    composition,diversity,lengths=[],[],[]
    for stat in rows(ROOT/"results/generation_statistics.tsv"):
        target,seed=stat["target"],int(stat["seed"])
        seqs=[r["sequence"] for r in rows(ROOT/f"data/generated/{target}/markov_{seed}.tsv.gz")]
        bases=Counter("".join(seqs));n=sum(bases.values())
        for b in "ACGT":composition.append(dict(target=target,method="markov",seed=seed,base=b,nucleotide_count=bases[b],fraction=bases[b]/n))
        for length,count in sorted(Counter(map(len,seqs)).items()):lengths.append(dict(target=target,method="markov",seed=seed,length=length,candidates=count))
        rng=random.Random(cfg["master_seed"]+seed)
        values=[]
        for _ in range(min(1000,len(seqs)*(len(seqs)-1)//2)):
            i,j=rng.sample(range(len(seqs)),2)
            values.append(Levenshtein.normalized_similarity(seqs[i],seqs[j]))
        diversity.append(dict(target=target,method="markov",seed=seed,sampled_pairs=len(values),sampling_seed=cfg["master_seed"]+seed,
             mean_pair_edit_identity=sum(values)/len(values),algorithm="1000_uniform_index_pairs_with_replacement_no_self_pairs; descriptive_QC_not_primary_endpoint"))
    table(ROOT/"results/nucleotide_composition.tsv",composition)
    table(ROOT/"results/generation_diversity.tsv",diversity)
    table(ROOT/"results/generation_lengths.tsv",lengths)


if __name__=="__main__":
    main()
