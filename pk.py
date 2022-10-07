#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
for paired-kmers

requires Python version >= 3.10 and NumPy

by Wubin Qu <quwubin@gmail.com>
2022-10-06 19:04:37
"""

import argparse
import numpy as np
import sys
import logging
import random

# ref: https://edinburgh-genome-foundry.github.io/easy_dna/_modules/easy_dna/random_sequences.html
def random_dna_sequence(length, gc_share=None, probas=None, seed=None):
    """Return a random DNA sequence ("ATGGCGT...") with the specified length.

    Parameters
    ----------

    length
      Length of the DNA sequence.

    gc_share
      The GC content of the random sequence, as a fraction (for example,
      0.3 for 30%). Overwrites `probas`.

    probas
      Frequencies for the different nucleotides, for instance
      ``probas={"A":0.2, "T":0.3, "G":0.3, "C":0.2}``.
      If not specified, all nucleotides are equiprobable (p=0.25).

    seed
      The seed to feed to the random number generator. When a seed is provided
      the random results depend deterministically on the seed, thus enabling
      reproducibility.
    """
    if seed is not None:
        np.random.seed(seed)
    if gc_share is not None:
        g_or_c = gc_share / 2.0
        not_g_or_c = (1 - gc_share) / 2.0
        probas = {"G": g_or_c, "C": g_or_c, "A": not_g_or_c, "T": not_g_or_c}
    if probas is None:
        sequence = np.random.choice(list("ATCG"), length)
    else:
        bases, probas = zip(*probas.items())
        sequence = np.random.choice(bases, length, p=probas)
    return "".join(sequence)


def gen_test_data(args):
    """
    
    parser_gen_test_data = subparsers.add_parser('gen-test-data', help=u'gene test data')
    parser_gen_test_data.add_argument('-left', type=str, default="cttaaatagggacctgtatgaatggctc", help='left conserved sequence')
    parser_gen_test_data.add_argument('-right', type=str, default="actaccaaacctgcattaaaaatttcgg", help='right conserved sequence')
    parser_gen_test_data.add_argument('-mins', type=int, default=30, help='[default 30] min space')
    parser_gen_test_data.add_argument('-maxs', type=int, default=200, help='[default 200] max space')
    parser_gen_test_data.add_argument('-minl', type=int, default=29000, help='[default 29000] min sequence size')
    parser_gen_test_data.add_argument('-maxl', type=int, default=31000, help='[default 31000] max sequence size')
    parser_gen_test_data.add_argument('-column', type=int, default=100000, help='[default 100,000] column size')
    parser_gen_test_data.add_argument('-number', type=int, default=1, help='[default 1] sequence number to generate')
    parser_gen_test_data.add_argument('-o', type=str, default="out.fa", help='output file name')
    parser_gen_test_data.set_defaults(func=gen_test_data)
    """
    seq_size_stop = args.maxs-len(args.left)-len(args.right)
    with open(args.o, 'w') as fh:
        for i in range(1, args.number+1):
            insert_size = random.randint(args.mins, seq_size_stop)
            init_seq_size = random.randint(args.minl, args.maxl)
            
            column_list = []
            if init_seq_size <= args.column:
                column_list.append(init_seq_size)
            else:
                column_list = [args.column] * (((init_seq_size - 1) // args.column) + 1 )

            seq_name = f"random_seq_{i}"
            seq = ""
            
            for seq_size in column_list:
                insert_pos = random.randint(1, seq_size - insert_size)
                left_part_size = insert_pos
                right_part_size = seq_size - (insert_pos + insert_size)
                seq += random_dna_sequence(left_part_size) + args.left + random_dna_sequence(insert_size) + args.right + random_dna_sequence(right_part_size)
                
            fh.write(f">{seq_name}\n{seq}\n")
    
    
def main():
    """Main."""
    parser = argparse.ArgumentParser(prog='./pk.py', description='commands for paired-kmers')

    subparsers = parser.add_subparsers(help='random data')
    subparsers.required = True

    parser_gen_test_data = subparsers.add_parser('gen-test-data', help=u'generate test data')
    parser_gen_test_data.add_argument('-left', type=str, default="cttaaatagggacctgtatgaatggctc", help='left conserved sequence')
    parser_gen_test_data.add_argument('-right', type=str, default="actaccaaacctgcattaaaaatttcgg", help='right conserved sequence')
    parser_gen_test_data.add_argument('-mins', type=int, default=30, help='[default 30] min space')
    parser_gen_test_data.add_argument('-maxs', type=int, default=200, help='[default 200] max space')
    parser_gen_test_data.add_argument('-minl', type=int, default=29000, help='[default 29,000] min sequence size')
    parser_gen_test_data.add_argument('-maxl', type=int, default=31000, help='[default 31,000] max sequence size')
    parser_gen_test_data.add_argument('-column', type=int, default=100000, help='[default 100,000] column size')
    parser_gen_test_data.add_argument('-number', type=int, default=1, help='[default 1] sequence number to generate')
    parser_gen_test_data.add_argument('-o', type=str, default="out.fa", help='output file name')
    parser_gen_test_data.set_defaults(func=gen_test_data)

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    logging.basicConfig(level=logging.DEBUG,
                        filename=f"run.log",
                        format='%(asctime)s %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S',
                        filemode='a')

    # define a Handler which writes INFO messages or higher to the sys.stderr
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)

    logging.info(f"{' '.join(sys.argv)}")

    args.func(args)


if __name__ == '__main__':
    main()
