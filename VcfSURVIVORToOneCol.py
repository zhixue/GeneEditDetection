#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/8/1 9:16 AM
    @Usage: python3 VcfSURVIVORToOneCol.py method_merged_survivor.vcf > new.vcf
"""
import sys
# need resort #ID=Chr, records
only_sort_flag = 0

with open(sys.argv[1]) as f:
    anno_lines = []
    temp_lines = []
    record_lines = dict()
    for line in f:
        if line.startswith('#'):
            if 'ID=Chr10' in line or \
                'ID=Chr11' in line or \
                'ID=Chr12' in line:
                temp_lines += [line]
                continue
            elif line.startswith('##ALT='):
                if temp_lines != []:
                    anno_lines += temp_lines
                    temp_lines = []
                anno_lines += [line]
            elif line.startswith('#CHROM'):
                temp = line.rstrip().split('\t')
                if len(temp[9:]) > 3:
                    only_sort_flag = 1
                    anno_lines += ['\t'.join(temp) + '\n']
                    continue
                temp = temp[:10]
                anno_lines += ['\t'.join(temp) + '\n']
            else:
                anno_lines += [line]
        else:
            temp = line.rstrip().split('\t')
            chrn = temp[0]
            if not chrn in record_lines:
                record_lines[chrn] = []
            if only_sort_flag:
                record_lines[chrn] += [line]
            else:
                gt_col = ''
                support_method = ''
                for i in range(9,len(temp)):
                    if not temp[i].startswith('.'):
                        if gt_col == '':
                            gt_col = temp[i]
                        if i == 9:
                            support_method += 'Delly,'
                        elif i == 10:
                            support_method += 'Manta,'
                        elif i == 11:
                            support_method += 'Smoove'
                info = temp[7].rstrip(';') + ';SUPP_METHOD=' + support_method.rstrip(',')
                temp[7] = info
                newline = '\t'.join(temp[:9] + [gt_col])
                record_lines[chrn] += [newline]

for an_l in anno_lines:
    print(an_l.rstrip())

for chi in range(1,13):
    chrn = 'Chr' + str(chi)
    if chrn not in record_lines:
        continue
    for rec_l in record_lines[chrn]:
        print(rec_l.rstrip())




