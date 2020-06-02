from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
from collections import Counter, defaultdict
import csv
from itertools import product
import math
import os


def letter_to_number(seq):
    conversion = {'A': '0', 'C': '1', 'G': '2', 'T': '3'}
    return ''.join(conversion[l] for l in seq)


def main(args):
    # Load each reference and strip newlines
    alphabet = '0123' if args.int else 'ACTG'
    source = open(args.source, 'r')
    references = []
    ref = ''
    emit = False if args.leave_out else True
    for line in source:
        if line.startswith('>'):
            if ref:
                if emit:
                    references.append(ref)
                    if args.reverse:
                        reverse_ref  = ref[::-1]
                        references.append(reverse_ref)
                    emit = False if args.leave_out else True
                ref = ''
            if args.leave_out and args.leave_out not in line:
                emit = True
            continue
        ref += line.strip()
    references.append(ref)
    if args.reverse:
        reverse_ref = ref[::-1]
        references.append(reverse_ref)
    source.close()

    # Initialize all words with 1 count
    word_length = args.word_length
    counter = defaultdict(Counter)
    #for word in (''.join(x) for x in product(alphabet, repeat=word_length)):
    #    counter[word[:-1]][word[-1]] = 1

    # Count word frequencies
    for ref in references:
        for i in range(len(ref) - word_length):
            prefix = ref[i:i+word_length-1]
            end_char = ref[i+word_length]
            if args.int:
                try:
                    prefix = letter_to_number(prefix)
                    end_char = letter_to_number(end_char)
                except KeyError:
                    continue
            counter[prefix][end_char] += 1

    hapax = 0
    for prefix in counter.values():
        for c in prefix.values():
            if c == 1:
                hapax += 1

    print(f'Hapax legomenom percentage: {100 * (hapax / (len(counter) * 4))}')

    # Calculate log probabilities
    probs = {}
    """ cache for sum of counts of words with prefix. Ex:
    >>> counter[AA]=2, counter[AC]=1, counter[AG]=1, counter[AT]=1
    >>> prefix_sum_cache[A] = 5
    """
    for prefix in counter.keys():
        prefix_sum = sum(counter[prefix].values())
        for end_char in counter[prefix].keys():
            prob = counter[prefix][end_char] / prefix_sum

            # skip .25 probabilities. In most cases these are OOV words
            if prob == 0.25:
                continue
            if args.log:
                prob = math.log(prob)
            probs[prefix+end_char] = prob

    # Write log probabilities to file destination
    file_name = '_'.join((args.name, str(args.word_length), 'gram'))
    if args.int:
        file_name += '_int'
    if args.log:
        file_name += '_log'
    if args.leave_out:
        file_name += f'_WO{args.leave_out}'
    file_name += '.lm'
    with open(os.path.join(args.dest, file_name), 'w+', newline='') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerows(probs.items())


def argparser():
    parser = ArgumentParser(
        formatter_class=ArgumentDefaultsHelpFormatter,
        add_help=False
    )
    parser.add_argument("source")
    parser.add_argument("dest")
    parser.add_argument('--word_length', type=int, default=12)
    parser.add_argument('--int', action='store_true')
    parser.add_argument('--log', action='store_true')
    parser.add_argument('--name', type=str, default='')
    parser.add_argument('--leave_out', type=str)
    parser.add_argument('--reverse', action='store_true')

    return parser
