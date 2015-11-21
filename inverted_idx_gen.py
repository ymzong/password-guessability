# Inverse-Index Builder for Markov Model Probabilities
#
# Input: Markov stat generated from statgen.py
# Output: Map from probability level to prefixes and transitions
#
# For 08-731 F15
# Authors: Derek Tzeng (dtzeng), Yiming Zong (yzong)

import json                             # For JSON I/O
import os                               # For tilde path
import math                             # For log, round
from collections import defaultdict     # Easier index mgmt

# Input statgen file prefix
STATGEN_PREFIX = os.path.expanduser("~/password-guessability/data/markov_chain_")
# Prefix for inverted index output
INDEX_PREFIX = "prob_index_"

def calc_scaling(stat, prob_levels):
    """
    Calculate the scaling factors C1 and C2
    Goal: ln(c1 * prob_min + c2) = 0; ln(c1 * prob_max + c2) = -prob_levels
    """
    print("Run One: Calculating scaling factors...")
    Transition = stat[0]
    StartCount = stat[1]

    # First obtain min probability and max probability
    prefix_count = sum(StartCount.itervalues())
    min_p = min(StartCount.itervalues()) * 1.0 / prefix_count
    max_p = max(StartCount.itervalues()) * 1.0 / prefix_count

    for prefix in Transition:
        prefix_count = sum(Transition[prefix].itervalues())
        min_p = min(min_p, min(Transition[prefix].itervalues()) * 1.0 / prefix_count)
        max_p = max(max_p, max(Transition[prefix].itervalues()) * 1.0 / prefix_count)

    # Solve the system of equations for c1 and c2
    c1 = (1 - math.exp(-prob_levels)) / (max_p - min_p)
    c2 = 1 - c1 * max_p
    print("Scaling parameters: c1={}, c2={}\n".format(c1, c2))

    return (c1, c2)


def build_invidx(statgen_f, prob_levels):
    with open(statgen_f, 'r') as f:
        stat = json.load(f)
    (c1, c2) = calc_scaling(stat, prob_levels)

    print("Run Two: Building inverted index for prefixes\n")
    Transition = stat[0]
    StartCount = stat[1]
    
    PrefixIdx = defaultdict(lambda: [])
    prefix_count = sum(StartCount.itervalues())
    for prefix in StartCount:
        scaled_p = math.log(c1 * StartCount[prefix] * 1.0 / prefix_count + c2)
        PrefixIdx[-int(round(scaled_p))].append(prefix)

    print("Run Three: Building inverted index for transitions\n")
    TransIdx = defaultdict(lambda: defaultdict(lambda: []))
    for prefix in Transition:
        prefix_count = sum(Transition[prefix].itervalues())
        for next_chr in Transition[prefix]:
            p = Transition[prefix][next_chr] * 1.0 / prefix_count
            scaled_p = math.log(c1 * p + c2)
            TransIdx[prefix][-int(round(scaled_p))].append(next_chr)

    return (TransIdx, PrefixIdx)


if __name__ == "__main__":
    prob_levels = int(raw_input("Probability levels (<= 10): "))
    k = int(raw_input("Order of Markov Model (<= 5): "))
    raw_input("Pleaes make sure the inverted index is `" + STATGEN_PREFIX + str(k) + "`...")
    print("Building inverted index...")
    Idx = build_invidx(STATGEN_PREFIX + str(k), prob_levels)

    OUTPUT_FILE = INDEX_PREFIX + str(k) + "_" + str(prob_levels)
    print("Writing output to {}...".format(OUTPUT_FILE))
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(Idx, f, sort_keys=True, indent=4)
    print("Done!")

