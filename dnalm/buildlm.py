from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import Counter
import csv
import math
import numpy as np
import os


def main(args):
    # Load each reference and strip newlines
    source = open(args.source, 'r')
    references = []
    ref = ''
    for line in source:
        if line.startswith('>'):
            if ref:
                references.append(ref)
                ref = ''
            continue
        ref += line.strip()
    references.append(ref)
    source.close()

    # Count word frequencies
    word_length = args.word_length
    counter = Counter()
    for ref in references:
        counter.update(Counter([ref[i:i+word_length] for i in range(len(ref) - word_length)]))

    print(f'Hapax legomenom percentage: {100 * (len([c for c in counter.values() if c == 1])) / len(counter)}')

    # Calculate log probabilities
    log_probs = {}
    prefix_sum_cache = {}
    """ cache for sum of counts of words with prefix. Ex:
    >>> counter[AA]=2, counter[AC]=1, counter[AG]=1, counter[AT]=1
    >>> prefix_sum_cache[A] = 5
    """
    filtered_count = 0
    for word in counter.keys():
        # Filter word with length above 8
        #if len(word) > 8:
        #    # These words needs at least 3 entries in counter
        #    if counter[word] <= 3:
        #        filtered_count += 1
        #        continue
        word_prefix = word[:-1]
        if prefix_sum_cache.get(word_prefix):
            prefix_sum = prefix_sum_cache[word_prefix]
        else:
            same_prefix_words = [w for w in counter.keys() if len(w) == len(word) and w.startswith(word_prefix)]
            prefix_sum = sum([counter[w] for w in same_prefix_words])
            prefix_sum_cache[word_prefix] = prefix_sum

        log_probs[word] = math.log(counter[word] / prefix_sum)

    # Write log probabilities to file destination
    with open(os.path.join(args.dest, f'{args.word_length}_gram.txt'), 'w+', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(log_probs.items())


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("source")
    parser.add_argument("dest")
    parser.add_argument("--word_length", type=int, default=12)

    return parser
