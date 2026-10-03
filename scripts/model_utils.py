from itertools import product
import math
import numpy as np
from common import settings


def enrichment(c_late, total_late, c_early, total_early, vocabulary, alpha=0.5):
    if min(total_late, total_early, vocabulary) <= 0 or alpha <= 0:
        raise ValueError("Positive depths, vocabulary and pseudocount required")
    return math.log2(((c_late+alpha)/(total_late+alpha*vocabulary)) /
                     ((c_early+alpha)/(total_early+alpha*vocabulary)))


def kmer_features(sequences, sizes=(1, 2, 3)):
    vocabulary = ["".join(p) for k in sizes for p in product("ACGT", repeat=k)]
    result = []
    for s in sequences:
        result.append([sum(s[i:i+len(k)] == k for i in range(len(s)-len(k)+1)) / max(1, len(s)-len(k)+1) for k in vocabulary])
    return np.asarray(result, dtype=float)


def assay_strength(value, endpoint):
    if endpoint == "KD":
        if value <= 0:
            raise ValueError("KD must be positive")
        return -math.log10(value)
    if endpoint == "SPR_end_of_injection_response":
        return value
    raise ValueError("Endpoint direction must be declared")


def binary_label(value):
    if str(value) not in ["0", "1"]:
        raise ValueError("Only published 0/1 labels accepted")
    return int(value)
