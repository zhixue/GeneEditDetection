#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/6/21 2:20 PM
    @Usage: python3 RebuildRef.py ref.fa ins.vcf TDNA.fa > new.fa
"""

import sys

def modify_fa(rawfa, ins_pos_dict, TDNA_dict):
    poss = [x for x in ins_pos_dict.keys()]
    rawfa_split = []
    if len(poss) == 1:
        rawfa_split = [rawfa[:poss[0]],
                       rawfa[poss[0]:]]
    else:
        last_pos = 0
        for i in range(len(poss)):
            pos = poss[i]
            if i == 0:
                rawfa_split += [rawfa[:pos]]
                last_pos = pos
            elif i == len(poss) - 1:
                rawfa_split += [rawfa[last_pos:pos],
                                rawfa[pos:]]
                last_pos = pos
            else:
                rawfa_split += [rawfa[last_pos:pos]]
                last_pos = pos

    tdnafa_split = []
    for i in range(len(poss)):
        pos = poss[i]
        tdnafa_split += [TDNA_dict[ins_pos_dict[pos]]]
    newseq = ''
    for i in range(len(poss)):
        newseq += rawfa_split[i].upper() + tdnafa_split[i].lower()
        if i == len(poss) - 1:
            newseq += rawfa_split[i+1].upper()
    return newseq



tdna_dict = dict()
with open(sys.argv[3]) as f:
    for line in f:
        if line.startswith('>'):
            seqid = line.rstrip().split(' ')[0][1:]
            tdna_dict[seqid] = ''
        else:
            tdna_dict[seqid] += line.rstrip()

ins_dict = dict()
with open(sys.argv[2]) as f:
    for line in f:
        if line.startswith('#'):
            continue
        temp = line.rstrip().split('\t')
        chrn = temp[0]
        pos = int(temp[1])
        seq = temp[2].split('_')[0]
        if not chrn in ins_dict:
            ins_dict[chrn] = dict()
        ins_dict[chrn][pos] = seq


seqid = ''
with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('>'):
            if seqid != '':
                if seqid in ins_dict:
                    #try:
                    newseqfa = modify_fa(seqfa, ins_dict[seqid], tdna_dict)
                    #except:
                    #    print(ins_dict[seqid])
                    #    exit()
                    print('>' + seqid + ' modify loc: ' + str(ins_dict[seqid]))
                else:
                    print('>' + seqid)
                    newseqfa = seqfa
                print(newseqfa)
            seqid = line.rstrip().split(' ')[0][1:]
            seqfa = ''
        else:
            seqfa += line.rstrip()
    # last one
    if seqid in ins_dict:
    #    try:
        newseqfa = modify_fa(seqfa, ins_dict[seqid], tdna_dict)
    #    except:
    #        print(ins_dict[seqid])
    #        exit()
    else:
        newseqfa = seqfa
    print('>' + seqid + ' modify loc: ' + str(ins_dict[seqid]))
    print(newseqfa)

