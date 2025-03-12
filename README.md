# TDNA Insertion infer

---
### Introduction
Infer TDNA insertion with known sites and known TDNA sequence

### Requirement
* samtools
* pysam (in python3)

### Usage
```shell
usage: INSBreakpointInfer.py [-h] -b <region.bed> -r <reference.given.region.sam> -i <TDNA.sam> -o <output_directory>

Infer the breakpoint and filter the reads, with given regions from reference genomes and TDNA sequences, checking the mapping position and CIGAR

options:
  -h, --help            show this help message and exit
  -b <region.bed>, --bed_path <region.bed>
                        Path of give region
  -r <reference.given.region.sam>, --reference_sam_path <reference.given.region.sam>
                        Path of sam of given region in reference
  -i <TDNA.sam>, --tdna_sam_path <TDNA.sam>
                        Path of sam in TDNA
  -o <output_directory>, --output_dir <output_directory>
                        Path of output
```
### Quick start with example data
Here we show an example of 100x read with 1%, 5% and 15% events TDNA insertion (tdna) in **rice** genome (WT:5000-15000).
* Step 1: Preparing TDNA insertion region bed with TDNA name
```shell
cd Example
mkdir test
cd test
echo -e "WT\t5000\t15000\ttdna" > test.bed
```

* Step 2: Preparing alignment results with regions

The bed has four coloums: **ChrID, start, end, TDNAname**

The region (end-start) in bed should **not** be too large, 10000 bp is recommanded.
```shell
for prefix in 1 5 15;
do
    refbam=../100x/100x_${prefix}_ref.bam
    tdnabam=../100x/100x_${prefix}_tdna.bam

  samtools index $refbam
  samtools view -h -M -L test.bed -O SAM -o 100x_${prefix}_ref_region.sam $refbam

  samtools index $tdnabam
  samtools view -h -F 4 -O SAM -o 100x_${prefix}_tdna.sam $tdnabam
done
```

* Step 3: Running INSBreakpointInfer.py

The results are in dir **100x_1_output**, **100x_5_output** and **100x_15_output**
```shell
for prefix in 1 5 15;
do
  python3 ../../INSBreakpointInfer.py -b test.bed -r 100x_${prefix}_ref_region.sam -i 100x_${prefix}_tdna.sam -o 100x_${prefix}_output
done
```
