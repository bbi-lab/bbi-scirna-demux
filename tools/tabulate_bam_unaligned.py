#!/usr/bin/env python3

import sys

# Sentinel-P5F02-P7F07_12|Sentinel|F02|F07|P08-H10_LIG184|TGCGAAAC        4       *       0       0       *       *       0       0       AGGCAGAGGCAGGCGGCTTTCTGAGTTCAAGGCCAGCCGGGTCTACAAAGTGAGNTCCANGCCCGCCAGGCCTATACAGAGAAA    CCCCCCCCCC-CCCCCCCCCCCCCCCCCCCCCCCCC;CCCCCCCCC;C-CCC;C#CCCC#CCC-CCC;CCCCCC--C-C;CCCC    sS:Z:AAGTTTGAAAGTGAAAACAATAAAAATGTGCGAAAC       sQ:Z:CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


#
# Program: tabulate_bam.py
# Purpose: make a table of read information from a rna_rtlig_demux BAM file.
# Notes:
#   o  this program is meant for testing by comparing the output to the
#      output from tabulate_fastq.py.
#   o  the output files must be sorted
#
for line in sys.stdin:
  line = line.rstrip()
  parts1 = line.split('\t')
  name = parts1[0]
  seq = parts1[9]
  qual = parts1[10]
  parts2 = name.split('|')
  p7 = parts2[2]
  p5 = parts2[3]
  rtlig = parts2[4]
  umi = parts2[5]
  print('%s|%s|%s|%s|%s|%s' % (p7, p5, rtlig, umi, seq, qual))
