# Higher-order Markov Model Stat Generator
#
# Input: Sanitized password file with count:
#        - first line is the (0-)index of delimiter
#        - following lines contain "${count} ${word}"
#   The input file is read from ../data/input/dataset-ascii.csv
#
# Output: The resulting starting, transition, and ending probabilities
# are written to local files in json format, with file name:
#   ../data/probs/#{value of k}_#{smoothing technique}_#{probability type}.json
#
# For 08-731 F15
# Authors: Derek Tzeng (dtzeng), Yiming Zong (yzong)

import json                             # For JSON I/O
from collections import defaultdict     # For easier count mgmt
import os                               # For path expansion
import string                           # For string constants
import itertools                        # For getting (k-1) grams

# Current directory of script
CURRENT_DIR = os.path.dirname(os.path.realpath('__file__'))
# Input password file
PASSWD_FILE = os.path.join(CURRENT_DIR, "../data/input/dataset-ascii.csv")
# Number of passwords between updating to stdout
UPDATE_INTERVAL = 1200000
# Valid chars in password
ALPHABET = string.digits + string.ascii_letters
ALPHABET_SIZE = len(ALPHABET)

# Smoothing constant
SMOOTH_DELTA = 0.01

# Range for k-gram size
K_GRAM_RANGE = xrange(2,  # lowest k
                      5 + 1,  # highest k + 1
                      )


def file_len(fname):
    """
    Given a text file, return the number of lines it has.
    """

    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def train_markov(StartCount, Transition, EndCount, k, passwd, freq):
    """
    Train markov model with a password (k-gram mode)
    """

    # Sanity check: if password too short / out of alphabet, ignore
    l = len(passwd)
    if l < k - 1:
        return
    for c in passwd:
        if c not in ALPHABET:
            return

    StartCount[passwd[:k - 1]] += freq          # Record prefix

    EndCount[passwd[-(k - 1):]] += freq         # Record termination

    # Regular transition cases
    for i in xrange(l - (k - 1)):
        substr = passwd[i:i + k - 1]
        next_chr = passwd[i + k - 1]
        Transition[substr][next_chr] += freq

    return


def build_markov_count(fname, k):
    """
    Train a higher-order Markov model on k-grams based on the
    password set. Returns the Markov chain as a tuple of three
    dictionaries.
    First dictionary maps starting (k-1)-gram to count; Second
    dictionary maps previous (k-1) grams to next chars then to
    count, i.e: { "foobar" : { 'a' : 12, 'b' : 9, ... }, ... }
    """

    StartCount = defaultdict(lambda: 0)
    Transition = defaultdict(lambda: defaultdict(lambda: 0))
    EndCount = defaultdict(lambda: 0)

    row_count = 0
    total_row = file_len(fname) - 1
    with open(PASSWD_FILE, 'r') as f:
        delimiter = int(f.readline())
        for l in f:
            word = (l[delimiter:-1]).strip()
            freq = int(l[:delimiter])
            train_markov(StartCount, Transition, EndCount, k, word, freq)

            row_count += 1
            if row_count % UPDATE_INTERVAL == 0:
                print("Progress: {} out of {} processed".format(
                    row_count, total_row))

    return (StartCount, Transition, EndCount)


def report_probability(StartCount, MidCount, EndCount, is_smoothed, delta):
    smooth_str = "additive" if is_smoothed else "none"

    # For starting probability
    tmp_dict = defaultdict(lambda: 0)
    dict_sum = sum(StartCount.values())
    new_sum = dict_sum + delta * ALPHABET_SIZE ** (k - 1)
    for pref in StartCount:
        tmp_dict[pref] = (StartCount[pref] + delta) * 1.0 / new_sum
    if is_smoothed:
        tmp_dict[""] = delta * 1.0 / new_sum        # Pseudo-count for non-existent prefix

    outfile = os.path.join(
        CURRENT_DIR, "../data/probs/{}_{}_start.json".format(k, smooth_str))
    print("Writing output to {}...".format(outfile)),
    with open(outfile, 'w') as f:
        json.dump(tmp_dict, f, sort_keys=True, indent=4)
    print("Done!")
    del(tmp_dict)

    # For ending probability
    tmp_dict = defaultdict(lambda: 0)
    dict_sum = sum(EndCount.values())
    new_sum = dict_sum + delta * ALPHABET_SIZE ** (k - 1)
    for suff in EndCount:
        tmp_dict[suff] = (EndCount[suff] + delta) * 1.0 / new_sum
    if is_smoothed:
        tmp_dict[""] = delta * 1.0 / new_sum        # Pseudo-count for non-existent suffix

    outfile = os.path.join(
        CURRENT_DIR, "../data/probs/{}_{}_end.json".format(k, smooth_str))
    print("Writing output to {}...".format(outfile)),
    with open(outfile, 'w') as f:
        json.dump(tmp_dict, f, sort_keys=True, indent=4)
    print("Done!")
    del(tmp_dict)

    # For transition probability
    tmp_dict = defaultdict(lambda: defaultdict(lambda: 0))
    for pref in MidCount:
        dict_sum = sum(MidCount[pref].values())
        new_sum = dict_sum + delta * ALPHABET_SIZE
        for next_chr in MidCount[pref]:
            tmp_dict[pref][next_chr] = MidCount[pref][next_chr] * 1.0 / new_sum

        if is_smoothed:
            tmp_dict[pref][""] = delta * 1.0 / new_sum  # Ditto

    outfile = os.path.join(
        CURRENT_DIR, "../data/probs/{}_{}_mid.json".format(k, smooth_str))
    print("Writing output to {}...".format(outfile)),
    with open(outfile, 'w') as f:
        json.dump(tmp_dict, f, sort_keys=True, indent=4)
    print("Done!")
    del(tmp_dict)


def additive_smooth_ends(d, k, delta):
    """
    Smoothes a (begin / end) frequency dictionary in place by adding
    pseudocount delta to every entry.
    """
    grams = (''.join(i) for i in itertools.product(ALPHABET, repeat=k - 1))

    for s in grams:
        d[s] += delta


def additive_smooth_middle(d, delta):
    """
    Smoothes a (transition) frequency dictionary in place by adding
    pseudocount delta to every entry.
    """

    grams = (''.join(i) for i in itertools.product(ALPHABET, repeat=k - 1))

    for s in d:
        for next_char in ALPHABET:
            d[s][next_char] += delta

# Main Routine
if __name__ == "__main__":
    print("Warning: For k=5, more than 6GB of memory will be used!")
    print("         You don't want to go beyond k=5 unless you have specialized hardware!")
    for k in K_GRAM_RANGE:
        print(
            "Processing {}-grams and without smoothing...".format(k))
        (StartCount, MidCount, EndCount) = build_markov_count(PASSWD_FILE, k)
        report_probability(StartCount, MidCount, EndCount, False, 0)
        report_probability(StartCount, MidCount, EndCount, True, SMOOTH_DELTA)
