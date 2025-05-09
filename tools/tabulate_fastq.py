#!/usr/bin/env python3

#
# Program: tabulate_fastq.py
# Purpose: make a table of read information from a rna_rtlig_demux BAM file.
# Notes:
#   o  run bbi-dmux/bin/make_sample_fastqs.py on the fastq files that were
#      demultiplexed by bcl-convert in a bbi-scirna-demux run.
#   o  this program is meant for testing by comparing the output to the
#      output from tabulate_bam.py. 
#   o  the output files must be sorted
#
#   o  typical script for running bbi-dmux/bin/make_sample_fastqs.py, which
#      was lifted from a bbi-dmux work_demux/*/*/.command.sh script
#
## #!/bin/bash -ue
## set -ueo pipefail
## 
## mkdir demux_out
## 
## source /net/gs/vol1/home/bge/git/bbi-dmux/load_pypy_env_reqs.sh
## source /net/gs/vol1/home/bge/git/bbi-dmux/bin/pypy_env/bin/activate
## 
## pypy3 /net/gs/vol1/home/bge/git/bbi-dmux/bin/make_sample_fastqs.py --run_directory /net/bbi/vol1/data/regression_tests/sciRNAseq/data_sources/illumina_runs/240514_VH00979_217_AAFJH2MM5/         --read1 <(zcat ../001_067_014_S19_L001_R1_001.fastq.gz) --read2 <(zcat ../001_067_014_S19_L001_R2_001.fastq.gz)         --file_name 001_067_014_S19_L001_R1_001.fastq.gz --sample_layout good_sample_sheet.csv         --p5_cols_used 2 2 --p7_rows_used E F         --p5_wells_used 0 --p7_wells_used 0         --pcr_index_pair_file 0         --rt_barcode_file default         --p5_barcode_file /net/bbi/vol2/home/bge/bbi/tests/nobackup/RNA3-72-a/p5.bmartin.20220927.txt         --p7_barcode_file default         --lig_barcode_file default         --multi_exp "0"         --output_dir ./demux_out --level 3
## 
## deactivate
## 
## pigz -p 8 demux_out/*.fastq
## 

# @Sentinel-P5F02-P7F07_1|Sentinel|F02|F07|P08-H10_LIG184|TGCGAAAC
# AGGCAGAGGCAGGCGGCTTTCTGAGTTCAAGGCCAGCCGGGTCTACAAAGTGAGNTCCANGCCCGCCAGGCCTATACAGAGAAA
# +
# CCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC;C-CCC;C#CCCC#CCC-CCC;CCCCCC--C-C;CCCC
# @Sentinel-P5F02-P7F07_2|Sentinel|F02|F07|P08-H10_LIG32|ATCCTGCA
# CCCTTTATCATGACTTTGATTTTCTAAGAGTGAAAGTACAAAGCAATTACAAAGGTTTTCCAGTAACAAGCATGACTGATGTAA
# +
# CCCCCCCCCCCC;CCCCCCCCCCCCCCC;CCC-CCCCCC--;CCC;CCCC-C-C;CCCCCCCC;CCC-C;CCCCCCCCCC;;CC
# @Sentinel-P5F02-P7F07_3|Sentinel|F02|F07|P05-H11_LIG323|GTATCGAT
# GGAGAGAACTTGCAGGGCAGTGTGTATGTTGTAGACACTAAAGTGCAACATTCGCTCAANAAAGATGTATTGGTAGCTTTCTTT


import sys

for (i, line) in enumerate(sys.stdin):
  line = line.rstrip()
  if(i % 4 == 0):
    parts = line.split('|')
    p7 = parts[2]
    p5 = parts[3]
    rtlig = parts[4]
    umi = parts[5]
  elif(i % 4 == 1):
    seq = line
  elif(i % 4 == 3):
    qual = line
    print('%s|%s|%s|%s|%s|%s' % (p7, p5, rtlig, umi, seq, qual))
