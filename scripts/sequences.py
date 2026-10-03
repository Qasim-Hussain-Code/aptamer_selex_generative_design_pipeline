"""Edit-distance families include insertions and deletions in short variable regions."""
import hashlib
import math
from rapidfuzz.distance import Levenshtein


def normalize(sequence):
    return sequence.upper().replace("U", "T").replace(" ", "")


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def identity(a, b):
    return 1.0 - Levenshtein.distance(a, b) / max(len(a), len(b))


def families(sequences, threshold):
    seqs = sorted(set(sequences))
    parent = list(range(len(seqs)))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    for i, a in enumerate(seqs):
        for j in range(i):
            b = seqs[j]
            cutoff = math.floor((1 - threshold) * max(len(a), len(b)) + 1e-9)
            if abs(len(a) - len(b)) > cutoff:
                continue
            if Levenshtein.distance(a, b, score_cutoff=cutoff) <= cutoff:
                x, y = find(i), find(j)
                if x != y:
                    parent[max(x, y)] = min(x, y)
    components = {}
    for i, seq in enumerate(seqs):
        components.setdefault(find(i), []).append(seq)
    ids = {k: "fam_" + hashlib.sha256("\n".join(v).encode()).hexdigest()[:16] for k, v in components.items()}
    return {seq: ids[find(i)] for i, seq in enumerate(seqs)}


def deterministic_subset(sequences, size, seed):
    return sorted(set(sequences), key=lambda s: (hashlib.sha256(f"{seed}:{s}".encode()).digest(), s))[:size]


def assert_no_leakage(train, test, assignment):
    overlap = {assignment[s] for s in train} & {assignment[s] for s in test}
    if set(train) & set(test) or overlap:
        raise AssertionError("Held-out family appears in training")
