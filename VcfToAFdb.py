#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/7/31 3:04 PM
    @Usage: python3 VcfToAFdb.py RiceVarMap2.vcf > RiceVarMap2.AF.tsv
"""
import sys

var_long = dict()
var_long_id = dict()
last_pos_long = -1
with open(sys.argv[1]) as f:
    print('\t'.join([str(x) for x in ['#CHROM', "POS", "ID", "REF", "ALT", "REF_AF", "ALT_AF", "Miss_AF"]]))
    for line in f:
        if line.startswith('#') or line.startswith('CHR'):
            continue
        temp = line.rstrip().split()
        chrn = temp[0].replace('chr', 'Chr').replace('Chr0', 'Chr')
        try:
            x = int(chrn)
            chrn = 'Chr' + chrn
        except:
            pass
        pos = int(temp[1])
        # SNP-SEEK SNP
        if len(temp) == 2:
            vid = '.'
            ref = '.'
            alt_str = '.'

            r_af = '.'
            v_af = '.'
            miss_af = '.'
            print('\t'.join(
                [str(x) for x in [chrn, pos, vid, ref, alt_str, r_af, v_af, miss_af]]))
        # SNP-SEEK SV and INDEL
        elif len(temp) == 6:
            if not chrn in var_long:
                var_long[chrn] = dict()
                var_long_id[chrn] = dict()
                last_pos_long = -1
            if not pos in var_long[chrn]:
                var_long[chrn][pos] = []
                var_long_id[chrn][pos] = []
            var_long[chrn][pos] += temp[4]
            # INS without ';'
            if ';' in temp[3]:
                vid = temp[3].split(';')[0] + temp[5] + '_' + 'LEN' + temp[3].split(';')[1]
            else:
                vid = 'INS' + temp[5] + '_' + 'LEN' + str(max([len(x) for x in temp[3].split(',')]))
            if not vid in var_long_id[chrn][pos]:
                var_long_id[chrn][pos] += [vid]

        # VCF, RiceVarMap2
        else:

            vid = temp[2]
            ref = temp[3]
            alt_str = temp[4]
            alt = temp[4].split(',')

            gt_n = len(temp[9:])
            ref_gt_n = len([x for x in temp[9:] if x.startswith('0|0') or x.startswith('0/0')])
            unknown_gt_n = len([x for x in temp[9:] if x.startswith('.')])
            r_af = ref_gt_n / gt_n
            v_af = (gt_n - ref_gt_n - unknown_gt_n) / gt_n
            miss_af = unknown_gt_n / gt_n
            print('\t'.join([str(x) for x in [chrn, pos, vid, ref, alt_str, round(r_af, 4), round(v_af, 4), round(miss_af, 4)]]))

# SNP-SEEK SV and INDEL
gt_n = 3023
if var_long:
    for chrn in var_long:
        for pos in var_long[chrn]:
            r_af = '.'
            miss_af = '.'
            vid = ','.join(var_long_id[chrn][pos])
            ref = '.'
            alt_str = '.'
            v_af = len(var_long[chrn][pos]) / gt_n
            print('\t'.join(
                [str(x) for x in [chrn, pos, vid, ref, alt_str, r_af, round(v_af, 4), miss_af]]))
