#!/usr/bin/env python3

#
# 'scrunch' a well-per-row samplesheet to a
# sample-per-row samplesheet.
#
# Reads SampleSheet.csv and writes to stdout.
#
# Notes:
#   o  this program is intended to convert a well-oriented samplesheet
#      file to a sample-oriented file. That is, one line per well becomes
#      one line per sample.
#   o  however, it works only when the rows for each sample are in blocks.
#      This should be fixed for consistency.
#

import io
import sys
import re
import argparse

# RT Barcode,Sample ID,Reference Genome
# P5-A01,24.0284,Mouse
# P7-A02,24.0284,Mouse
# P6-A01,24.0284,Mouse
# P5-A02,24.0284,Mouse
# P7-A01,24.0284,Mouse
# P6-A02,24.0284,Mouse
# 
# P6-A03,24.0292,Mouse
# P5-A04,24.0292,Mouse
# P5-A03,24.0292,Mouse
# P6-A04,24.0292,Mouse
# P7-A03,24.0292,Mouse
# P7-A04,24.0292,Mouse

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make bcl-convert samplesheet file.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input CSV samplesheet filename (required string).')
  parser.add_argument('-o', '--output', required=True, default=None, help='Output CSV samplesheet filename (required string).')
  parser.add_argument('-7', '--p7_wells', required=False, default=None, help='P7 well range string.')
  parser.add_argument('-5', '--p5_wells', required=False, default=None, help='P5 well range string.')
  parser.add_argument('-r', '--p7_rows', required=False, default=None, help='P7 row range string.')
  parser.add_argument('-c', '--p5_columns', required=False, default=None, help='P5 column range string.')
  parser.add_argument('-v', '--version', required=False, default=None, help='Write version string to stdout.')
  args = parser.parse_args()

  # print('Fix me: wrong headers for p7_rows and p5_columns')

  pcr_wells = False
  if(args.p7_wells != None and args.p5_wells != None):
    p7_spec = args.p7_wells
    p5_spec = args.p5_wells
    pcr_wells = True
  elif(args.p7_rows != None and args.p5_columns != None):
    p7_spec = args.p7_rows
    p5_spec = args.p5_columns
    pcr_wells = False
  elif(args.p7_wells != None and args.p5_wells == None):
    p7_spec = args.p7_wells 
    p5_spec = 'none'
    pcr_wells = True
  elif(args.p7_wells == None and args.p5_wells != None):
    p7_spec = 'none'
    p5_spec = args.p5_wells 
    pcr_wells = True
  else:
    p7_spec = 'none'
    p5_spec = 'none'
    pcr_wells = False

  ofh = open(args.output, 'w')

  if(pcr_wells == True):
    print('rt_wells,sample_name,genome,p7_wells,p5_wells', file=ofh)
  else:
    print('rt_wells,sample_name,genome,p7_rows,p5_columns', file=ofh)

  old_sample_name = ''
  with open(args.input, 'r') as fh:
    for line in fh:
      if(re.match('RT Barcode,Sample ID,Reference Genome', line)):
        continue
      line = line.strip()
      if(len(line) == 0):
        continue
      toks = line.split(',')
      sample_name = toks[1]
      if(sample_name != old_sample_name):
        if(old_sample_name != ''):
          print('\"%s\",%s,%s,\"%s\",\"%s\"' % (well_string, old_sample_name, genome, p7_spec, p5_spec), file=ofh)
        well_string = ''
        old_sample_name = sample_name
        genome=toks[2]
      if(len(well_string) > 0):
        well_string += ','
      well_string += toks[0]
    print('\"%s\",%s,%s,\"%s\",\"%s\"' % (well_string, old_sample_name, genome, p7_spec, p5_spec), file=ofh)

