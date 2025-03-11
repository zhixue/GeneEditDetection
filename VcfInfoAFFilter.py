#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/8/5 3:08 PM
    @Usage: python3 VcfInfoAFFilter.py raw.vcf KeyAF . > out.line
"""
import sys

def string2dict(long_string, sep=';', eq='=', rm_quote=False):
    if rm_quote:
        long_string = long_string.replace('"', '').replace("'", '')
    long_string = long_string.replace('; ', ';')
    out_dict = dict()
    tmp = long_string.rstrip(sep).split(sep)
    for i in tmp:
        if len(i.split(eq)) == 2:
            key, value = i.split(eq)
        else:
            key = i.split(eq)[0]
            value = eq.join(i.split(eq)[1:])
        out_dict[key] = value
    return out_dict

infokey = sys.argv[2]
cutoff = sys.argv[3]
try:
    cutoff = float(cutoff)
except:
    cutoff = '.'

header = ''
with open(sys.argv[1]) as fin:
    for line in fin:
        if line.startswith('#'):
            header += line
        else:
            temp = line.rstrip().split('\t')
            info = string2dict(temp[7])
            if infokey in info:
                if cutoff == '.':
                    print(line.rstrip)
                elif float(info[infokey]) <= cutoff:
                    print(line.rstrip)
                else:
                    continue