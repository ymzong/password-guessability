# Checkpoint builder based on discrete probabilities
#
# Input: Discrete probabilities from ../data/levels/*_*_*.json.
#
# Output: For each `interval` number of attempts, output the current guess
#         given the length and level.
#   The outputs are written to ../data/checkpoints/${k}_${smoothing}/${len}_${level}.out
#
# For 08-731 F15
# Authors: Derek Tzeng (dtzeng), Yiming Zong (yzong)

import json                             # For JSON I/O
import os                               # For path expansion
import sys                              # For argv and exit
import math                             # For log, round
from collections import defaultdict     # Easier index mgmt
import string                           # Use string constants
import itertools                        # Fancy list functions

# Current directory of script
CURRENT_DIR = os.path.dirname(os.path.realpath('__file__'))
# Input directory for level files
INPUT_PREFIX = os.path.join(CURRENT_DIR, "../data/levels/")
# Output directory for checkpoint file
CHECKPOINT_PREFIX = os.path.join(CURRENT_DIR, "../data/checkpoints/")

# Valid chars in password
ALPHABET = string.digits + string.ascii_letters
ALPHABET_SIZE = len(ALPHABET)

# Maximal level
MAX_LEVEL = 10

# How often do we want to checkpoint?
UPDATE_FREQUENCY = 10000

# Guess count so far
GUESS_COUNT = 0

# Scaling factors for probability (curated for k=3 with smoothing)
c1 = 1.11439835558
c2 = 4.53959983946e-05
# Default lvl for next char
NEXT_CHR_LVL = -int(round(math.log(c1 / ALPHABET_SIZE + c2)))

# List of passwords to be included in the checkpoint
CHECKPOINT = []


def enumerate_passwords(k, smoothing, l, total_level, freq):
    if l < k - 1:
        print("ERROR: Length of password too short (< k-1)!")
        return  # Nothing to do here :)

    ##################
    # Load input files
    file_name = INPUT_PREFIX + "{}_{}_start.json".format(k, smoothing)
    with open(file_name, 'r') as f:
        global start_lvl
        start_lvl = json.load(f)
    start_lvl = {int(key): val for key, val in start_lvl.iteritems()}
    global start_tokens
    start_tokens = set(itertools.chain.from_iterable(
        start_lvl.values()))   # Record all starting (k-1)-grams

    file_name = INPUT_PREFIX + "{}_{}_end.json".format(k, smoothing)
    with open(file_name, 'r') as f:
        global end_lvl
        end_lvl = json.load(f)
    end_lvl = {int(key): val for key, val in end_lvl.iteritems()}
    global end_tokens
    end_tokens = set(itertools.chain.from_iterable(end_lvl.values())
                     )       # Record all ending (k-1)-grams

    file_name = INPUT_PREFIX + "{}_{}_mid.json".format(k, smoothing)
    with open(file_name, 'r') as f:
        global mid_lvl
        mid_lvl = json.load(f)
    for prefix in mid_lvl:
        mid_lvl[prefix] = {int(key): val for key, val in mid_lvl[prefix].iteritems()}
    global mid_tokens
    mid_tokens = {key: set(itertools.chain.from_iterable(val.values()))
                  for key, val in mid_lvl.iteritems()}  # Record chars after (k-1)-gram

    ##################
    # Let's enumerate!
    dfs_passwords(l,            # total length,
                  k,            # k-gram model
                  0,            # next char index to enumerate
                  "",           # password so far
                  total_level,  # total level remaining
                  )


def dfs_passwords(l, k, next_idx, passwd, remaining_lvl):
    global GUESS_COUNT
    global CHECKPOINT

    # Finishing case
    if next_idx == l:   # End of password
        if remaining_lvl == 0:
            GUESS_COUNT += 1
            if GUESS_COUNT % UPDATE_FREQUENCY == 0:
                CHECKPOINT.append(passwd)
        return
    # Trim impossible cases
    elif MAX_LEVEL * (l - next_idx) < remaining_lvl:
        return
    # Last character
    elif next_idx == l - 1:
        suffix = passwd[-(k - 1):]
        if suffix not in mid_lvl:
            if remaining_lvl != NEXT_CHR_LVL:
                return
            for c in ALPHABET:
                GUESS_COUNT += 1
                if GUESS_COUNT % UPDATE_FREQUENCY == 0:
                    CHECKPOINT.append(passwd + c)
        elif remaining_lvl not in mid_lvl[suffix]:
            return  # Ehh, bad case
        for next_chr in mid_lvl[suffix][remaining_lvl]:
            if next_chr != "":
                GUESS_COUNT += 1
                if GUESS_COUNT % UPDATE_FREQUENCY == 0:
                    CHECKPOINT.append(passwd + next_chr)
            else:
                for c in ALPHABET:
                    if c in mid_tokens[suffix]:
                        continue
                    GUESS_COUNT += 1
                    if GUESS_COUNT % UPDATE_FREQUENCY == 0:
                        CHECKPOINT.append(passwd + c)
    # Initial case
    elif next_idx == 0:
        for init_level in xrange(0, min(remaining_lvl, MAX_LEVEL) + 1):
            if init_level not in start_lvl:
                continue
            for init_sequence in start_lvl[init_level]:
                if init_sequence != "":
                    # Normal case
                    dfs_passwords(l, k, k - 1, init_sequence, remaining_lvl - init_level)
                else:
                    # Wildcard case...
                    for init_seq in ("".join(k) for k in itertools.product(ALPHABET, repeat=k - 1)):
                        if init_seq in start_tokens:
                            continue
                        dfs_passwords(l, k, k - 1, init_seq, remaining_lvl - init_level)
    # Intermediate case
    else:
        prefix = passwd[-(k - 1):]
        if prefix not in mid_lvl:
            # Special case when we apply uniform probability to everything
            for c in ALPHABET:
                dfs_passwords(l, k, next_idx + 1, passwd + c, remaining_lvl - NEXT_CHR_LVL)
            return
        # Regular case
        for next_level in xrange(0, min(remaining_lvl, MAX_LEVEL) + 1):
            if next_level not in mid_lvl[prefix]:
                continue
            for next_chr in mid_lvl[prefix][next_level]:
                # Wildcard case...
                if next_chr == "":
                    for c in ALPHABET:
                        if c in mid_tokens[prefix]:
                            continue
                        dfs_passwords(l, k, next_idx + 1, passwd + c, remaining_lvl - next_level)
                # Normal case
                else:
                    dfs_passwords(l, k, next_idx + 1, passwd + next_chr, remaining_lvl - next_level)

if __name__ == "__main__":
    # Input handling
    try:
        K = int(sys.argv[1])
        SMOOTHING = sys.argv[2]
        IS_SMOOTHING = (SMOOTHING != 'none')
        LEN = int(sys.argv[3])
        TOTAL_LEVEL = int(sys.argv[4])
    except Exception:
        print("usage: checkpoint.py K smoothing_mode length total_level\n")
        sys.exit(1)

    OUTPUT_FILE = os.path.join(CURRENT_DIR,
                               "../data/checkpoints/{}_{}/{}_{}.out".format(
                                   K, SMOOTHING, LEN, TOTAL_LEVEL)
                               )
    print("This script will enumerate passwords with length {} and level {},".format(LEN, TOTAL_LEVEL))
    print("Output will be written to {}".format(OUTPUT_FILE))
    print("Enumerating......")

    enumerate_passwords(K, SMOOTHING, LEN, TOTAL_LEVEL, UPDATE_FREQUENCY)

    # Write output to checkpoint file
    with open(OUTPUT_FILE, 'w') as f:
        for passwd in CHECKPOINT:
            f.write("{}\n".format(passwd))
        f.write("\n")
        f.write(str(GUESS_COUNT) + "\n")

    print("Checkpointing finished! Total passwords: {}.".format(GUESS_COUNT))
