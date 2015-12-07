# Inverted-Index Builder for discretized Markov Probabilities
#
# Input: Markov stat generated from statgen.py, located in:
#   ../data/probs/*.json
#
# Output:
#   - For starting / ending probabilities, return the inverted index
#     from a probability level to prefixes and suffixes;
#   - For transition probabilities, return the mapping from each (k-1)-
#     gram to an inner mapping from level to next characters.
#   The outputs are written to ../data/levels/*_*_*.json.
#
# For 08-731 F15
# Authors: Derek Tzeng (dtzeng), Yiming Zong (yzong)

import json                             # For JSON I/O
import os                               # For path expansion
import math                             # For log, round
from collections import defaultdict     # Easier index mgmt

# Current directory of script
CURRENT_DIR = os.path.dirname(os.path.realpath('__file__'))
# Input directory for probability files
INPUT_PREFIX = os.path.join(CURRENT_DIR, "../data/probs/")
# Output directory for level files
OUTPUT_PREFIX = os.path.join(CURRENT_DIR, "../data/levels/")

# Range for k-gram size
K_GRAM_RANGE = xrange(2, 5 + 1)
# Range for smoothing options
SMOOTHING_LIST = ("none",       # No smoothing applied
                  "additive",   # Additive smoothing
                  )

# Max-level to use for discretization
LEVEL = 10


def calc_scaling(stat, two_layer=False):
    """
    Calculate the scaling factors C1 and C2
    Goal: ln(c1 * prob_min + c2) = 0; ln(c1 * prob_max + c2) = -prob_levels
    """
    print("Calculating scaling factors...")

    # First obtain min probability and max probability
    (min_p, max_p) = (2, 0)
    if two_layer:
        for prefix in stat:
            min_p = min(min_p, min(stat[prefix].itervalues()))
            max_p = max(max_p, max(stat[prefix].itervalues()))
    else:
        min_p = min(stat.itervalues())
        max_p = max(stat.itervalues())

    # Solve the system of equations for c1 and c2
    c1 = (1 - math.exp(-LEVEL)) / (max_p - min_p)
    c2 = 1 - c1 * max_p
    print("Scaling parameters: c1={}, c2={}".format(c1, c2))
    return (c1, c2)


def report_levels(k, s):
    """
    Build the discrete index based on input file of probabilities.
    """
    ##################################
    # Build level index for StartProbs
    file_name = INPUT_PREFIX + "{}_{}_start.json".format(k, s)
    with open(file_name, 'r') as f:
        start_p = json.load(f)
    (c1, c2) = calc_scaling(start_p)

    start_lvl = defaultdict(lambda: [])
    for prefix in start_p:
        scaled_p = -int(round(math.log(c1 * start_p[prefix] + c2)))
        start_lvl[scaled_p].append(prefix)

    outfile = OUTPUT_PREFIX + "{}_{}_start.json".format(k, s)
    print("Writing output to {}...".format(outfile))
    with open(outfile, 'w') as f:
        json.dump(start_lvl, f, sort_keys=True, indent=4)
    # Cleanup memory!
    del(start_p)
    del(start_lvl)

    ################################
    # Build level index for EndProbs
    file_name = INPUT_PREFIX + "{}_{}_end.json".format(k, s)
    with open(file_name, 'r') as f:
        end_p = json.load(f)
    (c1, c2) = calc_scaling(end_p)

    end_lvl = defaultdict(lambda: [])
    for suffix in end_p:
        scaled_p = -int(round(math.log(c1 * end_p[suffix] + c2)))
        end_lvl[scaled_p].append(suffix)

    outfile = OUTPUT_PREFIX + "{}_{}_end.json".format(k, s)
    print("Writing output to {}...".format(outfile))
    with open(outfile, 'w') as f:
        json.dump(end_lvl, f, sort_keys=True, indent=4)
    # Cleanup memory!
    del(end_p)
    del(end_lvl)

    ################################
    # Build level index for MidProbs
    file_name = INPUT_PREFIX + "{}_{}_mid.json".format(k, s)
    with open(file_name, 'r') as f:
        mid_p = json.load(f)
    (c1, c2) = calc_scaling(mid_p, True)    # Signal double-layer dictionary

    mid_lvl = defaultdict(lambda: defaultdict(lambda: []))
    for prefix in mid_p:
        for next_chr in mid_p[prefix]:
            scaled_p = -int(round(math.log(c1 * mid_p[prefix][next_chr] + c2)))
            mid_lvl[prefix][scaled_p].append(next_chr)

    outfile = OUTPUT_PREFIX + "{}_{}_mid.json".format(k, s)
    print("Writing output to {}...".format(outfile))
    with open(outfile, 'w') as f:
        json.dump(mid_lvl, f, sort_keys=True, indent=4)
    # Cleanup memory!
    del(mid_p)
    del(mid_lvl)


if __name__ == "__main__":
    print("Building inverted index with max-lvl {} from input ../data/probs/*".format(LEVEL))
    for k in K_GRAM_RANGE:
        for s in SMOOTHING_LIST:
            print("Processing {}-gram probabilities with smoothing option {}...".format(k, s))
            report_levels(k, s)
    print("Done!")
