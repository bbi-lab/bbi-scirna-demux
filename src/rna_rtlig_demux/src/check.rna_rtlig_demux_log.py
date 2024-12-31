#!/usr/bin/env python3

import io
import sys
import re
import argparse
import json


#
# check that the whitelist and mismatch sequences are distinct
#
def check_for_distinct_sequences(rt_counter_vec):
  print('check that whitelist and mismatches sequences are distinct')
  sequence_dict = {}
  number_sequences_checked = 0
  for i, barcode_counter in enumerate(rt_counter_vec):
    sequence_count = sequence_dict.setdefault(barcode_counter['whitelist_sequence'], 0)
    sequence_count += 1
    sequence_dict[barcode_counter['whitelist_sequence']] = sequence_count
#    print('i: %d  barcode count: %d' % (i, sequence_dict[barcode_counter['whitelist_sequence']]))
    number_sequences_checked += 1

    for j, mismatch_string in enumerate(barcode_counter['mismatch_sequences']):
      mobj = re.match(r'([ACGTN]+):([0-9]+)$', mismatch_string)
      mismatch_seq = mobj.group(1)
      sequence_count = sequence_dict.setdefault(mismatch_seq, 0)
      sequence_count += 1
      sequence_dict[mismatch_seq] = sequence_count
      number_sequences_checked += 1

  for sequence in sequence_dict.keys():
    sequence_count = sequence_dict[sequence]
    if(sequence_count != 1):
      print('oops: %s  %d' % (sequence, sequence_count))

  print('%d  sequences checked' % (number_sequences_checked))
  print('done')
  print()


#
# Check that there is exactly one mismatch in the
# mismatch sequences.
#
def check_mismatch_sequences(rt_counter_vec):
  print('check mismatch sequences have one mismatch')
  number_whitelist_sequences_checked = 0
  number_mismatch_sequences_checked = 0
  for i, barcode_counter in enumerate(rt_counter_vec):
    whitelist_sequence = barcode_counter['whitelist_sequence']
    mismatch_sequence_list = barcode_counter['mismatch_sequences']
    for mismatch_string in mismatch_sequence_list:
      mobj = re.match(r'([ACGTN]+):([0-9]+)$', mismatch_string)
      mismatch_seq = mobj.group(1)
      mismatch_read_count = int(mobj.group(2))
      num_mismatch = 0
      for i in range(len(whitelist_sequence)):
        if(whitelist_sequence[i] != mismatch_seq[i]):
          num_mismatch += 1
      if(num_mismatch != 1):
        print('yikes')
      number_mismatch_sequences_checked += 1
    number_whitelist_sequences_checked += 1
  print('  %d whitelist sequences checked' % (number_whitelist_sequences_checked))
  print('  %d mismatch sequences checked' % (number_mismatch_sequences_checked))
  print('done')
  print()


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


def index_to_well( well_index, across_row_first ):
  if( well_index < 0 ):
    return( (0, 'none'))
  nrow = 8
  ncol = 12
  ipl = int( well_index / 96 )
  i96 = well_index - ipl * 96
  if across_row_first:
      well_row = chr(65 + int(i96 / ncol))
      well_col = (i96 % ncol) + 1
  else:
      well_row = chr(65 + (i96 % nrow))
      well_col = int(i96 / nrow) + 1

#    well_id = 'P%d-%s%s' % (ipl + 1, well_row, pad_well_col(well_col, zero_pad_col, id_length))
  well_id = '%s%s' % ( well_row, pad_well_col( well_col, True, 2 ) )
  return( (ipl, well_id ) )
    

#
# Compare well and sample names in barcode_counter to
# those in the original (one row per well) samplesheet
# file.
#
def check_rt_wells_and_sample_names(rt_counter_vec, filename):
  print('check that barcode counter sample names match samplesheet well sample names')
  index_well_list = []
  index_well_list.append('Undetermined')
  samplesheet_dict = {}
  for well_index in range(8*96):
    ipl, well_id = index_to_well( well_index, True)
    well_name = 'P%02d-%s' % (ipl + 1, well_id)
    samplesheet_dict.setdefault(well_name, 'Undetermined')
    index_well_list.append(well_name)

  fh = open(filename, 'r')
  for row_string in fh:
    row_toks = row_string.strip().split(',')
    if(row_toks[0] == 'RT Barcode'):
      continue
    well_name = row_toks[0]
    mobj = re.match(r'P([0-9]+)-([A-H][0-9]+)', well_name)
    well_name = 'P%02d-%s' % (int(mobj.group(1)), mobj.group(2))
    sample_name = row_toks[1]
    samplesheet_dict[well_name] = sample_name

  num_wells_checked = 0
  for i, barcode_counter in enumerate(rt_counter_vec):
    whitelist_sequence = barcode_counter['whitelist_sequence']
    well_name = barcode_counter['well_name']
    sample_name = barcode_counter['sample_name']
    if(sample_name != samplesheet_dict[well_name]):
      print('wrong sample name: %s  %s  %s' % (sample_name, well_name, samplesheet_dict[well_name]))
    if(barcode_counter['well_name'] != index_well_list[int(barcode_counter['well_index'])]):
      print('wrong well index: %d  %s' % (well_index, barcode_counter['well_name'])) 
    num_wells_checked += 1
  print('  %d well-sample name pairs checked' % (num_wells_checked))
  print('done')
  print()


def check_rt_wells_and_barcode_seqs(rt_counter_vec, filename):
  print('check that barcode counter whitelist sequences match barcodes in rt.txt')
  barcode_dict = {}
  fh = open(filename, 'r')
  for row_string in fh:
    row_toks = row_string.strip().split()
    well_name = row_toks[0]
    barcode_sequence = row_toks[1]
    barcode_dict.setdefault(well_name, barcode_sequence)

  num_barcodes_checked = 0
  for i, barcode_counter in enumerate(rt_counter_vec):
    whitelist_sequence = barcode_counter['whitelist_sequence']
    well_name = barcode_counter['well_name']
    if(whitelist_sequence != barcode_dict[well_name]):
      print('wrong barcode sequence (%s) for well %s' % (whitelist_sequence, well_name))
    num_barcodes_checked += 1
  print(' %d barcode sequences checked' % (num_barcodes_checked))
  print('done')
  print()


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make bcl-convert samplesheet file.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON rna_rtlig_demux log filename (required string).')
  args = parser.parse_args()

  fh = open(args.input, 'r')
  log_json = json.load(fh)
  rt_counter_vec = log_json['rt_counter_vec']
  check_for_distinct_sequences(rt_counter_vec)
  check_mismatch_sequences(rt_counter_vec)
  samplesheet_filename = '/home/brent/work/data_sets/RNA3-065-a/SampleSheet.csv'
  check_rt_wells_and_sample_names(rt_counter_vec, samplesheet_filename)
  rt_barcode_filename = '/home/brent/git/bbi-dmux/bin/barcode_files/rt.txt'
  check_rt_wells_and_barcode_seqs(rt_counter_vec, rt_barcode_filename)


