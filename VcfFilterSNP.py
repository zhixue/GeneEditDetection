#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/7/31 8:42 AM
    @Usage: python3 VcfFilterSNP.py -i -o -e
"""
import argparse
import sys
import os

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Get vcf of SNP''')
    parser.add_argument('-i', '--in_vcf', metavar='<in.vcf>', help='Path of in.vcf', type=str,
                        required=True)
    parser.add_argument('-o', '--out_vcf', metavar='<out.vcf>', help='Path of out.vcf', type=str,
                    required=True)
    parser.add_argument('-e', '--exclude_SNP', metavar='<out.vcf>', help='Exclude SNP (default:0)', type=str,
                    default=0)

    args = vars(parser.parse_args())
    if os.path.exists(args["out_vcf"]):
        print('# Warning: Output vcf exists!')
        #exit()
    fout = open(args["out_vcf"], 'w')
    with open(args["in_vcf"]) as fin:
        for line in fin:
            if line.startswith('#'):
                fout.write(line)
            else:
                temp = line.rstrip().split('\t')
                chrn = temp[0]
                pos = int(temp[1])
                ref = temp[3]
                alt = temp[4].split(',')
                if '<NON_REF>' in alt:
                    alt.remove('<NON_REF>')
                snp_flag = 0
                if len(ref) == 1 and max([len(x) for x in alt]) == 1:
                    snp_flag = 1
                if args["exclude_SNP"]:
                    if not snp_flag:
                        fout.write(line)
                else:
                    if snp_flag:
                        fout.write(line)


