#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    @Author: Hongzhang Xue
    @Modified: 2024/6/21 8:30 AM
    @Usage: python3 INSBreakpointInfer.py -b region.bed -r region.sorted.sam -i TDNA.sorted.sam -o xx_dir
    # xx_dir:
    # xx_dir.tsv xx_dir.breakpoint.sam
"""
import argparse
import sys
from collections import Counter
import statistics
import pysam
import os

# before: samtools view -h -L xx.bed -O SAM -o xx.sam xx.bam
## xx.bed
#   Chr Start End   TDNA_ID
Omit_read_freq = 1

def coord2region(region_dict, chrn, coord):
    if not chrn in region_dict:
        return None
    for region in region_dict[chrn]:
        if region[0] <= coord <= region[1]:
            return region
    return None

def cigarsplit(st):
    st_list = []
    current_unit = ''
    for i in range(len(st)):
        if st[i] in [str(x) for x in range(0,10)]:
            current_unit += st[i]
        else:
            current_unit += st[i]
            st_list += [current_unit]
            current_unit = ''
    return st_list

def bool_paired(reada, readb):
    if reada.query_name != readb.query_name:    # exclude different reads
        return 0
    if reada.flag & 64 != readb.flag & 64 : # not both left/ right read
        return 0
    if not (reada.cigartuples and readb.cigartuples):  # unmapped end in paired-end reads
        return 0
    if len(readb.cigartuples) == 1:  # exclude only include M
        return 1
    if reada.flag & 16 == readb.flag & 16:  # mapped to the same strand
        if reada.cigartuples[0][0] == 0 and readb.cigartuples[-1][0] == 0 or \
                (reada.cigartuples[-1][0] == 0 and readb.cigartuples[0][0] == 0):   # 32S43M & 32M43S
            return 1
        if len(reada.cigartuples) == 3 and reada.cigartuples[-1][0] != 0 and reada.cigartuples[-1][1] <= 3 and readb.cigartuples[0][0] == 0 or \
                (reada.cigartuples[0][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[-1][0] != 0 and readb.cigartuples[-1][1] <= 3) or \
                (len(reada.cigartuples) == 3 and reada.cigartuples[0][0] != 0 and reada.cigartuples[0][1] <= 3 and readb.cigartuples[-1][0] == 0) or \
                (reada.cigartuples[-1][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[0][0] != 0 and readb.cigartuples[0][1] <= 3):  # 41S30M2S & 41MXXX | 2S30M41S & XXX41M
            return 1
    else:
        if reada.cigartuples[0][0] == 0 and readb.cigartuples[0][0] == 0 or \
                (reada.cigartuples[-1][0] == 0 and readb.cigartuples[-1][0] == 0):  # 32S43M & 43S32M
            return 1
        if len(reada.cigartuples) == 3 and reada.cigartuples[-1][0] != 0 and reada.cigartuples[-1][1] <= 3 and readb.cigartuples[-1][0] == 0 or \
                (reada.cigartuples[-1][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[-1][0] != 0 and readb.cigartuples[-1][1] <= 3) or \
                (len(reada.cigartuples) == 3 and reada.cigartuples[0][0] != 0 and reada.cigartuples[0][1] <= 3 and readb.cigartuples[0][0] == 0) or \
                (reada.cigartuples[0][0] == 0 and len(readb.cigartuples) == 3 and readb.cigartuples[0][0] != 0 and readb.cigartuples[0][1] <= 3):  # 41S30M2S & XXX41M | 2S30M41S & 41MXXX
            return 1

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='''Infer the breakpoint and filter the reads, 
    with given regions from reference genomes and TDNA sequences,
    checking the mapping position and CIGAR''')
    parser.add_argument('-b', '--bed_path', metavar='<region.bed>',
                        help='Path of give region', type=str, required=True)
    parser.add_argument('-r', '--reference_sam_path', metavar='<reference.given.region.sam>',
                        help='Path of sam of given region in reference', type=str, required=True)
    parser.add_argument('-i', '--tdna_sam_path', metavar='<TDNA.sam>',
                        help='Path of sam in TDNA', type=str, required=True)
    parser.add_argument('-o', '--output_dir', metavar='<output_directory>',
                        help='Path of output', type=str, required=True)
    args = vars(parser.parse_args())



    chr_region = dict()

    with open(args['bed_path']) as f:
        for line in f:
            temp = line.rstrip().split('\t')
            chrn = temp[0]
            start = int(temp[1])
            end = int(temp[2])
            if len(temp) > 3:
                tag = temp[3]
            else:
                tag = ''
            if not chrn in chr_region:
                chr_region[chrn] = dict()
            chr_region[chrn][(start, end)] = tag


    current_chr = ''
    current_region = (0, 0)
    left_points_coords = []
    left_points_reads = []
    right_points_coords = []
    right_points_reads = []

    outdir = args['output_dir']
    os.system('mkdir ' + outdir)
    fouttsv = open(outdir + '/' + outdir.split('/')[-1]  + '.tsv', 'w')
    foutsam = outdir + '/' + outdir.split('/')[-1]  + 'breakpoint.ref.sam'
    foutsam2 = outdir + '/' + outdir.split('/')[-1]  + 'breakpoint.TDNA.sam'

    fouttsv.write('\t'.join(['#Chr','(Start,End)', 'Tag', 'Breakpoint', 'InferredCoord',
                     'MedianCoord', 'STDCCoord', 'SupportReadTable']) + '\n')

    # load tdna sam
    read_tdna_dict = dict()
    read_tdna_pysam_dict = dict()
    #with open(args['tdna_sam_path']) as f:
    sam_tdna = pysam.AlignmentFile(args['tdna_sam_path'], 'r')
    for one in sam_tdna:
        readid = one.query_name
        flag = one.flag
        chrn = one.reference_name
        coord = one.reference_start
        cigar = one.cigarstring
        if not readid in read_tdna_dict:
            read_tdna_dict[readid] = []
        read_tdna_dict[readid] += [[flag, chrn, coord, cigar]]
        if not readid in read_tdna_pysam_dict:
            read_tdna_pysam_dict[readid] = []
        read_tdna_pysam_dict[readid] += [one]





    #with open(args['reference_sam_path']) as f:
    ref_filtered_read = []
    tdna_filtered_read = []

    sam_ref = pysam.AlignmentFile(args['reference_sam_path'], 'r')
    for one in sam_ref:
        readid = one.query_name
        flag = one.flag
        chrn = one.reference_name
        coord = one.reference_start
        cigar = one.cigarstring
        if cigar is None:
            continue
        if chrn.strip() == '':
            continue
        find_region = coord2region(chr_region, chrn, coord)
        if find_region is None:
            continue
        tdna_tag = chr_region[chrn][find_region].split('_')[0]
        #print(chrn,current_chr)
        #print(find_region)
        if current_chr!=chrn or find_region != current_region:
            if current_region != (0, 0):
                # summary
                left_points_coords_freq = Counter(left_points_coords)
                left_max_freq_coord = -1
                left_max_freq = 0
                left_sum_freq_Omit = 0
                left_sum_freq_0 = 0
                for point in left_points_coords_freq:
                    if left_points_coords_freq[point] > left_max_freq:
                         left_max_freq_coord = point
                         left_max_freq = left_points_coords_freq[point]
                    left_sum_freq_0 += left_points_coords_freq[point]
                    if left_points_coords_freq[point] > Omit_read_freq:
                        left_sum_freq_Omit += left_points_coords_freq[point]
                if len(left_points_coords) > 1:
                    left_points_coord_stdev = statistics.stdev(left_points_coords)
                    left_points_coord_median = statistics.median(left_points_coords)
                else:
                    if len(left_points_coords) > 0:
                        left_points_coord_median = left_points_coords[0]
                    else:
                        left_points_coord_median = 0
                    left_points_coord_stdev = 0

                right_points_coords_freq = Counter(right_points_coords)
                right_max_freq_coord = -1
                right_max_freq = 0
                right_sum_freq_Omit = 0
                right_sum_freq_0 = 0
                for point in right_points_coords_freq:
                    if right_points_coords_freq[point] > right_max_freq:
                        right_max_freq_coord = point
                        right_max_freq = right_points_coords_freq[point]
                    right_sum_freq_0 += right_points_coords_freq[point]
                    if right_points_coords_freq[point] > Omit_read_freq:
                        right_sum_freq_Omit += right_points_coords_freq[point]
                if len(right_points_coords) > 1:
                    right_points_coord_stdev = statistics.stdev(right_points_coords)
                    right_points_coord_median = statistics.median(right_points_coords)
                else:
                    if len(right_points_coords) > 0:
                        right_points_coord_median = right_points_coords[0]
                    else:
                        right_points_coord_median = 0
                    right_points_coord_stdev = 0

                fouttsv.write('\t'.join([str(x) for x in (current_chr, current_region, chr_region[current_chr][current_region], 'Left',
                                                      left_max_freq_coord,
                                                      left_points_coord_median, left_points_coord_stdev,
                                                      left_sum_freq_0,
                                                      left_sum_freq_Omit,
                                                      left_points_coords_freq,
                                                      )]) + '\n')
                fouttsv.write('\t'.join([str(x) for x in (current_chr, current_region, chr_region[current_chr][current_region], 'Right',
                                                      right_max_freq_coord,
                                                      right_points_coord_median, right_points_coord_stdev,
                                                      right_sum_freq_0,
                                                      right_sum_freq_Omit,
                                                      right_points_coords_freq,
                                                      )]) + '\n')
            current_chr = chrn
            current_region = find_region
            left_points_coords = []
            left_points_reads = []
            right_points_coords = []
            right_points_reads = []

        # match cigar
        #try:
        cigarlist = cigarsplit(cigar)
        pass_flag = 0
        if readid in read_tdna_dict:
            read_results = read_tdna_pysam_dict[readid]
            for res in read_results:
                # check flag, chrn, coord, cigar
                if tdna_tag == res.reference_name:
                    if bool_paired(one, res):
                        ref_filtered_read += [one]
                        tdna_filtered_read += [res]
                        #print(one, res)
                        pass_flag = 1
                    break
        #except:
        #    print(cigar)
        #    exit()
        # left breakpoint
        if pass_flag == 0:
            continue
        if cigarlist[-1].endswith('S'):
            left_move_bp = 0
            for k in range(len(cigarlist) -1):
                if cigarlist[k].endswith('M'):
                    left_move_bp += int(cigarlist[k][:-1])
                elif cigarlist[k].endswith('D'):
                    left_move_bp += int(cigarlist[k][:-1])
                elif cigarlist[k].endswith('I'):
                    left_move_bp -= int(cigarlist[k][:-1])
            left_points_coords += [coord + left_move_bp]
            left_points_reads += [readid]
            #foutsam.write(line)

        # right breakpoint
        if cigarlist[0].endswith('S'):
            right_move_bp = 0
            right_points_coords += [coord]
            right_points_reads += [readid]
            #foutsam.write(line)


    # finl one
    # summary
    left_points_coords_freq = Counter(left_points_coords)
    left_max_freq_coord = -1
    left_max_freq = 0
    left_sum_freq_Omit = 0
    left_sum_freq_0 = 0
    for point in left_points_coords_freq:
        if left_points_coords_freq[point] > left_max_freq:
            left_max_freq_coord = point
            left_max_freq = left_points_coords_freq[point]
        left_sum_freq_0 += left_points_coords_freq[point]
        if left_points_coords_freq[point] > Omit_read_freq:
            left_sum_freq_Omit += left_points_coords_freq[point]
    if len(left_points_coords) > 1:
        left_points_coord_stdev = statistics.stdev(left_points_coords)
        left_points_coord_median = statistics.median(left_points_coords)
    else:
        if len(left_points_coords) > 0:
            left_points_coord_median = left_points_coords[0]
        else:
            left_points_coord_median = 0
        left_points_coord_stdev = 0

    right_points_coords_freq = Counter(right_points_coords)
    right_max_freq_coord = -1
    right_max_freq = 0
    right_sum_freq_Omit = 0
    right_sum_freq_0 = 0
    for point in right_points_coords_freq:
        if right_points_coords_freq[point] > right_max_freq:
            right_max_freq_coord = point
            right_max_freq = right_points_coords_freq[point]
        right_sum_freq_0 += right_points_coords_freq[point]
        if right_points_coords_freq[point] > Omit_read_freq:
            right_sum_freq_Omit += right_points_coords_freq[point]
    if len(right_points_coords) > 1:
        right_points_coord_stdev = statistics.stdev(right_points_coords)
        right_points_coord_median = statistics.median(right_points_coords)
    else:
        if len(right_points_coords) > 0:
            right_points_coord_median = right_points_coords[0]
        else:
            right_points_coord_median = 0
        right_points_coord_stdev = 0

    fouttsv.write(
        '\t'.join([str(x) for x in (current_chr, current_region, chr_region[current_chr][current_region], 'Left',
                                    left_max_freq_coord,
                                    left_points_coord_median, left_points_coord_stdev,
                                    left_sum_freq_0,
                                    left_sum_freq_Omit,
                                    left_points_coords_freq,
                                    )]) + '\n')
    fouttsv.write(
        '\t'.join([str(x) for x in (current_chr, current_region, chr_region[current_chr][current_region], 'Right',
                                    right_max_freq_coord,
                                    right_points_coord_median, right_points_coord_stdev,
                                    right_sum_freq_0,
                                    right_sum_freq_Omit,
                                    right_points_coords_freq,
                                    )]) + '\n')

    # write sam
    #foutsam = outdir + '/' + outdir.split('/')[-1] + 'breakpoint.ref.sam'
    #foutsam2 = outdir + '/' + outdir.split('/')[-1] + 'breakpoint.TDNA.sam'
    out_ref_sam = pysam.AlignmentFile(foutsam, 'w', template=sam_ref)
    for oneR in ref_filtered_read:
        out_ref_sam.write(oneR)
    out_ref_sam.close()
    out_tdna_sam = pysam.AlignmentFile(foutsam2, 'w', template=sam_tdna)
    for oneT in tdna_filtered_read:
        out_tdna_sam.write(oneT)
    out_tdna_sam.close()
    sam_ref.close()
    sam_tdna.close()

