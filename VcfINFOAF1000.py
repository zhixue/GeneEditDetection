#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/7/31 10:23 AM
    @Usage: python3 VcfINFOAF.py raw.vcf RiceVarMap2.vcf SNP-Seek.SNP.txt > raw.addInfo.vcf
"""
# run VcfToAFdb first
# python3 VcfToAFdb.py RiceVarMap2.vcf > RiceVarMap2.AF.tsv

import sys
# SNP 0-base
# SV 1000-base
indel_max_dis = 1000

snp_seek_svsAF = dict()
with open(sys.argv[3]) as f:
    for line in f:
        if line.startswith('#'):
            continue
        temp = line.rstrip().split()
        chrn = temp[0]
        pos = int(temp[1])
        if not chrn in snp_seek_svsAF:
            snp_seek_svsAF[chrn] = dict()
        try:
            snp_seek_svsAF[chrn][pos] = float(temp[6])
        except:
            snp_seek_svsAF[chrn][pos] = '.'

ricevarmap2_snpsAF = dict()
ricevarmap2_indelsAF = dict()
with open(sys.argv[2]) as f:
    for line in f:
        if line.startswith('#'):
            continue
        temp = line.rstrip().split('\t')
        chrn = temp[0]
        pos = int(temp[1])
        ref = temp[3]
        alt = temp[4].split(',')
        snp_flag = 0
        indel_flag = 0
        if len(ref) == 1 and min([len(x) for x in alt]) == 1:
            snp_flag = 1
            if  max([len(x) for x in alt]) > 1:
                indel_flag = 1
        else:
            indel_flag = 1

        v_af = float(temp[6])
        if snp_flag:
            if not chrn in ricevarmap2_snpsAF:
                ricevarmap2_snpsAF[chrn] = dict()
            ricevarmap2_snpsAF[chrn][pos] = v_af
        if indel_flag:
            if not chrn in ricevarmap2_indelsAF:
                ricevarmap2_indelsAF[chrn] = dict()
            ricevarmap2_indelsAF[chrn][pos] = v_af


with open(sys.argv[1]) as f:
    for line in f:
        if line.startswith('#'):
            print(line.rstrip())
            continue
        temp = line.rstrip().split('\t')
        chrn = temp[0]
        pos = int(temp[1])
        ref = temp[3]
        alt = temp[4].split(',')
        snp_flag = 0
        indel_flag = 0
        if len(ref) == 1 and min([len(x) for x in alt]) == 1:
            snp_flag = 1
            if max([len(x) for x in alt]) > 1:
                indel_flag = 1
        else:
            indel_flag = 1

        info_str = temp[7]
        if snp_flag:
            if chrn in snp_seek_svsAF:
                if pos in snp_seek_svsAF[chrn]:
                    info_str = info_str.rstrip(';') + ';SNPSeekAF=' + str(snp_seek_svsAF[chrn][pos])
            if chrn in ricevarmap2_snpsAF:
                if pos in ricevarmap2_snpsAF[chrn]:
                    info_str = info_str.rstrip(';') + ';RiceVarMap2AF=' + str(round(ricevarmap2_snpsAF[chrn][pos], 4))
        if indel_flag:
            if chrn in ricevarmap2_indelsAF:
                min_dis = indel_max_dis
                min_dis_pos = -1
                for posRVM in ricevarmap2_indelsAF[chrn]:
                    if abs(pos - posRVM) <= indel_max_dis:
                        temp_dis = abs(pos - posRVM)
                        if temp_dis < min_dis:
                            min_dis_pos = posRVM
                    if posRVM - pos > indel_max_dis:
                        break
                if min_dis_pos != -1:
                    info_str = info_str.rstrip(';') + ';RiceVarMap2AF=' + str(
                        round(ricevarmap2_indelsAF[chrn][min_dis_pos], 4))

            if chrn in snp_seek_svsAF:
                min_dis = indel_max_dis
                min_dis_pos = -1
                for posSeek in snp_seek_svsAF[chrn]:
                    if abs(pos - posSeek) <= indel_max_dis:
                        temp_dis = abs(pos - posSeek)
                        if temp_dis < min_dis:
                            min_dis_pos = posSeek
                    if posSeek - pos > indel_max_dis:
                        break
                if min_dis_pos != -1:
                    info_str = info_str.rstrip(';') + ';SNPSeekAF=' + str(snp_seek_svsAF[chrn][min_dis_pos])


        temp[7] = info_str
        print('\t'.join(temp))





