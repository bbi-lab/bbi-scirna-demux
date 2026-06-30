#!/usr/bin/env python3

import sys

# 25.0464-P5none-P7B04_6039786|25.0464|none|B04|P05-A01_LIG4|ATGAACTA     16      chr1    21838341        255     54M     *       0       0       CTTGGGGCTGACTTATAGGTTTAGAGGTTCAGTCCATTATCATCAAGGTAGGAG  III9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9III   NH:i:1  HI:i:1  AS:i:53 nM:i:0

# 25.0464-P5none-P7P01-B04_12117817|25.0464|none|P01-B04|P05-A03_LIG57|TCGTCCCA  16       chr1    3111773 255     81M     *       0       0       AACATAGTACTTGAAGTAGTAGCCAGAGCAATTCGACAACAATAGGAGATCAAGGGGATATAAATTGGAAAAGAGGAAGAC       IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII99IIIIIIIIIIIIIIIIIIIIIIIIII      NH:i:1   HI:i:1  nM:i:0  AS:i:79 GX:Z:-  GN:Z:-  sM:i:0  CB:Z:AACGAATAAAATGCAAAACAAAAAAAAA       CY:Z:CCCCCCCCCCCCCCCCCCCCCCCCCCCC       UB:Z:TCGTCCCA   UY:Z:CCCCCCCC


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
  rnam = parts1[2]
  pos  = parts1[3]
  name = parts1[0]
  seq = parts1[9]
  qual = parts1[10]

  # edit name
  parts2 = name.split('|')
  parts3 = parts2[0].split('_')
  prefix = parts3[0]
  name_edited = prefix + '|' + parts2[1] + '|' + parts2[2] + '|' + parts2[3] + '|' + parts2[4] + '|' + parts2[5]
 
  print('%s %s %s %s %s' % (name_edited, seq, qual, rnam, pos))
