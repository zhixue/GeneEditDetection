#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/7/26 3:29 PM
    @Usage: python3 LongToCount.py raw.count sps.count > out.matrix
"""
import sys

sps_depth = dict()
with open(sys.argv[2]) as f:
    for line in f:
        temp = line.rstrip().split('\t')
        sps_depth[temp[0]] = float(temp[4])

event_depth = dict()
with open(sys.argv[1]) as f:
    for line in f:
        temp = line.rstrip().split('\t')
        if temp[3] not in event_depth:
            event_depth[temp[3]] = dict()
        event_depth[temp[3]][temp[0]] = float(temp[4]) / sps_depth[temp[0]]

print('Event\t' + '\t'.join(event_depth[temp[3]].keys()))
for evt in event_depth:
    s = evt
    for key in event_depth[evt]:
        s += '\t' + str(round(event_depth[evt][key],5))
    print(s)


