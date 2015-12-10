# Guess count calculator given discretized Markov probabilities and checkpoints
#
# Input:
#   - Discrete probabilities in ../data/levels/*_*_*.json.
#   - Password checkpoints in ../data/checkpoints/${k}_${smoothing}/${len}_${level}.out
#   - Input password
#
# Output: Guess number for the password, or BEYOND_THRESHOLD
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
import locale                           # For readable numeric output

locale.setlocale(locale.LC_ALL, '')

# Colorful shell output! :)


class color:
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    RED = '\033[91m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    END = '\033[0m'

# Valid chars in password
ALPHABET = string.digits + string.ascii_letters
ALPHABET_SIZE = len(ALPHABET)

# Checkpoint parameters
K = 3
SMOOTHING = "additive"
LVL_FACTOR = 2

# How often do we want to checkpoint?
CHECKPOINT_FREQUENCY = 10000

# Maximal level & length
MAX_LEVEL = 10
MAX_LENGTH = 12

# Scaling factors for probability (curated for k=3 with smoothing)
c1 = 1.11439835558
c2 = 4.53959983946e-05
# Default lvl for next char
NEXT_CHR_LVL = -int(round(math.log(c1 / ALPHABET_SIZE + c2)))

# Current directory of script
CURRENT_DIR = os.path.dirname(os.path.realpath('__file__'))
# Input directory for level files
LEVEL_PREFIX = os.path.join(CURRENT_DIR, "../data/levels/")
# Input directory for checkpoint files
CHECKPOINT_PREFIX = os.path.join(CURRENT_DIR, "../data/cps/{}_{}/".format(K, SMOOTHING))

# Available checkpoints
AVAILABLE_CP = {4: xrange(33), 5: xrange(35), 6: xrange(33), 7: xrange(
    28), 8: xrange(25), 9: xrange(23), 10: xrange(22), 11: xrange(21), 12: xrange(21)}


def load_levels(k, smoothing):
    ##################
    # Load input files
    global start_lvl
    global end_lvl
    global mid_lvl

    file_name = LEVEL_PREFIX + "{}_{}_start.json".format(k, smoothing)
    with open(file_name, 'r') as f:
        global start_lvl
        start_lvl = json.load(f)
    start_lvl = {int(key): val for key, val in start_lvl.iteritems()}
    global start_tokens
    start_tokens = set(itertools.chain.from_iterable(
        start_lvl.values()))   # Record all starting (k-1)-grams

    file_name = LEVEL_PREFIX + "{}_{}_end.json".format(k, smoothing)
    with open(file_name, 'r') as f:
        global end_lvl
        end_lvl = json.load(f)
    end_lvl = {int(key): val for key, val in end_lvl.iteritems()}
    global end_tokens
    end_tokens = set(itertools.chain.from_iterable(end_lvl.values())
                     )       # Record all ending (k-1)-grams

    file_name = LEVEL_PREFIX + "{}_{}_mid.json".format(k, smoothing)
    with open(file_name, 'r') as f:
        global mid_lvl
        mid_lvl = json.load(f)
    for prefix in mid_lvl:
        mid_lvl[prefix] = {int(key): val for key, val in mid_lvl[prefix].iteritems()}
    global mid_tokens
    mid_tokens = {key: set(itertools.chain.from_iterable(val.values()))
                  for key, val in mid_lvl.iteritems()}  # Record chars after (k-1)-gram


def decompose_password(pw, k):
    """
    Given a password, go through it char-by-char and return a tuple of tuples of the form
    (level, idx, prefix / char). level stands for the level of the prefix / next char;
    idx stands for the index of the option in the level index; prefix / char stand for the
    actual prefix / next char.

    The purpose of this is that we may compare two passwords with same length and total level
    with this tuple, in order to get their comparative order in the DFS sequence.
    """
    result = ()

    # Level for prefix
    prefix = pw[:(k - 1)]
    found = False
    index = -1
    for l in start_lvl:
        if prefix in start_lvl[l]:
            found = True
            index = start_lvl[l].index(prefix)
            break
    # Does not appear in index, must be wildcard
    if not found:
        for l in start_lvl:
            if "" in start_lvl[l]:
                index = start_lvl[l].index("")
                break
    assert(index >= 0)
    result += ((l, index, prefix),)

    # Level for intermediate characters
    for i in xrange(0, len(pw) - k + 1):
        prefix = pw[i: i + k - 1]
        # Prefix not in index -- all is wildcard case
        if prefix not in mid_lvl:
            result += ((NEXT_CHR_LVL, ALPHABET.index(pw[i + k - 1]), pw[i + k - 1]), )
            continue
        # Prefix in index -- find corresponding level
        found = False
        index = -1
        for l in mid_lvl[prefix]:
            if pw[i + k - 1] in mid_lvl[prefix][l]:
                found = True
                index = mid_lvl[prefix][l].index(pw[i + k - 1])
                break
        if not found:
            for l in mid_lvl[prefix]:
                if "" in mid_lvl[prefix][l]:
                    index = mid_lvl[prefix][l].index("")
                    break
        assert(index >= 0)
        result += ((l, index, pw[i + k - 1]), )
    return result


def skip_prev_cases(LEN, LVL):
    """
    Given (length, level) pair, go through and skip all previous cases with smaller or equal
    complexity value (len + lvl / beta). Sum up the number of skipped passwords along the way.
    """

    global guess_count

    max_complexity = LEN + LVL / LVL_FACTOR
    found = False
    for complexity in xrange(max_complexity + 1):
        for ln in xrange(min(AVAILABLE_CP.keys()), max(AVAILABLE_CP.keys()) + 1):
            for lvl in xrange((complexity - ln) * LVL_FACTOR, (complexity - ln + 1) * LVL_FACTOR):
                if lvl not in AVAILABLE_CP[ln]:
                    continue
                if (ln, lvl) == (LEN, LVL):
                    return True
                # Accumulate counts of prior cases
                # Note: For large files, only seek for the last line
                MAX_LINE = 100
                SEEK_END = 2
                f = open(CHECKPOINT_PREFIX + "{}_{}.out".format(ln, lvl))
                try:
                    f.seek(-MAX_LINE, SEEK_END)
                    guess_count += int(f.read(MAX_LINE).splitlines()[-1])
                except IOError:
                    guess_count += int(f.read().splitlines()[-1])
                f.close()
    return False


def binary_search(list, pw_components):
    """
    Binary search on the list of passwords to give the right-most one that gets
    enumerated no later than pw. This is done by comparing the "serialized
    passwords."
    """
    # Trivial case
    if len(list) == 0:
        return -1

    head = 0
    tail = len(list) - 1        # Inclusive
    while head <= tail:
        mid = (head + tail) / 2
        mid_pass = decompose_password(list[mid], K)
        if mid_pass > pw_components:
            tail = mid - 1
        else:
            head = mid + 1
    return head


def dfs_passwords(l, k, next_idx, passwd, remaining_lvl, lower_bound=None):
    """
    Nearly identical to the dfs function in checkpoint.py, except that we simplify
    it to discard the "one more char" case and include a lower_bound parameter for
    DFS pruning. Also, now the function returns True if the current run matches the
    desired password.
    """
    global guess_count
    global pw
    # Before anything else, validate the lower bound for previously chosen prefix / char.
    if lower_bound and next_idx > 0:
        prev_last_component = decompose_password(passwd, k)[-1]
        if lower_bound[0] < prev_last_component:
            lower_bound = None
        elif lower_bound[0] == prev_last_component:
            lower_bound = lower_bound[1:]
            if len(lower_bound) == 0:
                lower_bound = None
        else:   # prev_last_component < lower_bound
            return False    # Prune

    # Finishing case
    if next_idx == l:   # End of password
        if remaining_lvl == 0:
            guess_count += 1
            if passwd == pw:
                return True
        return False     # Bad case!
    # Trim impossible cases
    elif MAX_LEVEL * (l - next_idx) < remaining_lvl:
        return False
    # Initial case
    elif next_idx == 0:
        for init_level in xrange(0, min(remaining_lvl, MAX_LEVEL) + 1):  # !!!
            if init_level not in start_lvl:
                continue
            for init_sequence in start_lvl[init_level]:
                if init_sequence != "":
                    # Normal case
                    result = dfs_passwords(l, k, k - 1, init_sequence,
                                           remaining_lvl - init_level, lower_bound)
                    if result:
                        return True
                else:
                    # Wildcard case...
                    for init_seq in ("".join(k) for k in itertools.product(ALPHABET, repeat=k - 1)):
                        if init_seq in start_tokens:
                            continue
                        result = dfs_passwords(l, k, k - 1, init_seq,
                                               remaining_lvl - init_level, lower_bound)
                        if result:
                            return True
    # Intermediate case
    else:
        prefix = passwd[-(k - 1):]
        if prefix not in mid_lvl:
            # Special case when we apply uniform probability to everything
            for c in ALPHABET:
                result = dfs_passwords(l, k, next_idx + 1, passwd + c,
                                       remaining_lvl - NEXT_CHR_LVL, lower_bound)
                if result:
                    return True
            return False
        # Regular case
        for next_level in xrange(0, min(remaining_lvl, MAX_LEVEL) + 1):  # !!!
            if next_level not in mid_lvl[prefix]:
                continue
            for next_chr in mid_lvl[prefix][next_level]:
                # Wildcard case...
                if next_chr == "":
                    for c in ALPHABET:
                        if c in mid_tokens[prefix]:
                            continue
                        result = dfs_passwords(l, k, next_idx + 1, passwd +
                                               c, remaining_lvl - next_level, lower_bound)
                        if result:
                            return True
                # Normal case
                else:
                    result = dfs_passwords(l, k, next_idx + 1, passwd +
                                           next_chr, remaining_lvl - next_level, lower_bound)
                    if result:
                        return True
    return False    # No luck this time :(


if __name__ == "__main__":
    # Clean screen
    tmp = os.system("clear")

    # Load password level index
    load_levels(K, SMOOTHING)
    print("Parameters of model: k={}, smoothing={}".format(K, SMOOTHING))

    # Analyze password input
    pw = raw_input("Input password to guess -> " + color.UNDERLINE)
    print(color.END)
    for c in pw:
        if c not in ALPHABET:
            print("Only alpha-numeric passwords are supported in this version!\n")
            sys.exit(0)
    if len(pw) > MAX_LENGTH:
        print("Longest password supported is len = {}...\n".format(MAX_LENGTH))
        sys.exit(0)

    LEN = len(pw)
    components = decompose_password(pw, K)
    print("Password components: {}".format(components))
    LVL = sum(l for l, _, _ in components)
    print("Password length: {}; Total level: {}".format(LEN, LVL))

    # Skip over previous (length, level) cases
    guess_count = 0
    found = skip_prev_cases(LEN, LVL)
    if not found:
        print "Password is beyond our index space with {:n} passwords!\n".format(guess_count)
        sys.exit(0)

    # Go through current (length, level) case to find lower bound
    print "Skipped over {:n} passwords! :)\n".format(guess_count)
    with open(CHECKPOINT_PREFIX + "{}_{}.out".format(LEN, LVL)) as f:
        current_cp = f.read().splitlines()[:-2]     # Entire file excluding last two summary lines
    idx = binary_search(current_cp, components)
    guess_count += max(0, idx) * CHECKPOINT_FREQUENCY
    print("Using binary search to estimate guess count...")
    print "Narrowed search between {:n} and {:n}".format(guess_count, guess_count + CHECKPOINT_FREQUENCY)

    # Construct lower bound in DFS
    if idx == 0 or len(current_cp) == 0:
        lower_bound = None
    else:
        lower_bound = decompose_password(current_cp[idx - 1], K)

    dfs_passwords(LEN,            # total length,
                  K,            # k-gram model
                  0,            # next char index to enumerate
                  "",           # password so far
                  LVL,  # total level remaining
                  lower_bound   # for password searching
                  )
    print(color.BOLD + color.UNDERLINE + color.YELLOW +
          "\nGuess Count: {:n}\n".format(guess_count) + color.END * 3)
