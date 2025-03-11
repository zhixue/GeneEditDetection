#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/6/21 3:46 PM
    @Usage: python3 Fasta_SmallLeterToBed.py raw.fa > out.bed
"""

# lest split N: 1*N
# OUTPUT:
# chr , chr length, start N (1-base), end N (1-base,open), length of N

import sys

LEASTN = 1

ctg_name = ''
ctg_seq = ''

def n2bed(name,seq):
    # loc 0-base
    # [)
    # chr , chr length, start N (1-base), end N (1-base,open), length of N
    current_loc = 0

    end_loc = len(seq)
    if end_loc == 0:
        return
    while current_loc<=end_loc - 1:
        if not seq[current_loc] in ('a','c','g','t','n'):
            current_loc += 1
        else:
            n_end = current_loc + 1
            while n_end <= end_loc - 1:
                if seq[n_end] in ('a','c','g','t','n'):
                    n_end += 1
                    if n_end == end_loc:
                        nlength = n_end - current_loc
                        print('\t'.join([name.split()[0], str(current_loc), str(n_end), str(nlength)]))
                        return
                else:
                    nlength = n_end - current_loc
                    if nlength >= LEASTN:
                        print('\t'.join([name.split()[0], str(current_loc), str(n_end), str(nlength)]))
                    current_loc = n_end + 1

                    break





with open(sys.argv[1]) as fin:
    for line in fin:
        if line.startswith('>'):
            if ctg_name != '':
                # block compute
                n2bed(ctg_name,ctg_seq)

            ctg_name = line.rstrip()[1:]
            ctg_seq = ''
        else:
            ctg_seq += line.rstrip()
    # last one
    n2bed(ctg_name, ctg_seq)