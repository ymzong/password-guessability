# Higher-order Markov Model Stat Generator
# The result is written to a local file in json format
#
# For 08-731 F15
# Authors: Derek Tzeng (dtzeng), Yiming Zong (yzong)

import json                             # For JSON I/O
from collections import defaultdict     # For easier count mgmt
import os                               # For tilde path expansion


# Input password file
PASSWD_FILE = os.path.expanduser("~/password-guessability/data/rockyou.txt")
# Prefix for password stat output
STAT_FILE_PREFIX = "markov_chain_"
# Number of passwords between updating to stdout
UPDATE_INTERVAL = 200000
# Valid chars in password
ALPHABET = range(32, 127) 


def file_len(fname):
    """
    Given a text file, return the number of lines it has.
    """

    with open(fname) as f:
        for i, l in enumerate(f):
            pass
    return i + 1


def train_markov(Transition, StartCount, k, passwd):
    """
    Train markov model with a password (k-gram mode)
    """

    # Sanity check: if password too short / out of alphabet, ignore
    l = len(passwd)
    if l < k - 1:
        return
    for c in passwd:
        if ord(c) not in ALPHABET:
            return

    StartCount[passwd[0:k-1]] += 1          # Record prefix
    
    Transition[passwd[-(k-1):]][0] += 1     # Record termination

    # Regular transition cases
    for i in xrange(l - (k - 1)):
        substr = passwd[i:i + k - 1]
        next_chr = passwd[i + k - 1]
        Transition[substr][next_chr] += 1

    return

def build_markov(fname, k):
    """
    Train a higher-order Markov model on k-grams based on the
    password set. Returns the Markov chain as a tuple of two
    dictionaries.
    First dictionary is the transition stats: Key is the current
    state (a k-1 gram), and Value is another dictionary, mapping
    potential next char to number of occurences, i.e:
    { 'a' : 12, 'b' : 9, ... }.
    """

    Transition = defaultdict(lambda: defaultdict(lambda: 0))
    StartCount = defaultdict(lambda: 0)
    count = 0
    total_count = file_len(fname)
    print("Gathering statistics...")
    with open(fname) as f:
        for passwd in f:
            if count % UPDATE_INTERVAL == 0:
                print("Progress: {} out of {} processed".format(
                    count, total_count))
            count += 1
            train_markov(Transition, StartCount, k, passwd[:-1])
    return (Transition, StartCount) 

if __name__ == "__main__":
    k = int(raw_input("Order of Markov Model (<= 5): "))
    raw_input("Please make sure the password file is `"
        + PASSWD_FILE + "`...")
    
    print("Building Markov-model with {}-grams...".format(k))
    Stats = build_markov(PASSWD_FILE, k)

    OUTPUT_FILE = STAT_FILE_PREFIX + str(k)
    print("Writing output to {}...".format(OUTPUT_FILE))
    with open(OUTPUT_FILE, 'w') as f:
        json.dump(Stats, f, sort_keys=True, indent=4)
    print("Done!")

