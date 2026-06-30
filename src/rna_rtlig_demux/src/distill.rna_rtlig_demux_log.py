#!/usr/bin/env python3

import io
import sys
import re
import argparse
import json


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

  well_id = '%s%s' % ( well_row, pad_well_col( well_col, True, 2 ) )
  return( (ipl, well_id ) )
    

if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make bcl-convert samplesheet file.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON rna_rtlig_demux log filename (required string).')
  args = parser.parse_args()

  fh = open(args.input, 'r')
  log_json = json.load(fh)
  rt_counter_vec = log_json['rt_counter_vec']
  check_mismatch_sequences(rt_counter_vec)


