#!/usr/bin/env python3

import sys
import argparse


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


def index_to_well(well_index, across_row_first):
  if(well_index < 0):
    return((0, 'none'))
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


def base_convert(i, b):
  if(i == 0):
    return([0])
  result = []
  while i > 0:
    result.insert(0, i % b)
    i = i // b
  return result


def index_base_4_to_encoded_index(index_base_4):
  encoder_base = ['A', 'C', 'G', 'T']
  encoded_index_list = list()
  for i in index_base_4:
    encoded_index_list.append(encoder_base[i])
  encoded_index = ''.join(encoded_index_list).rjust(7, 'A')
  return(encoded_index)


if __name__ == '__main__':

  parser = argparse.ArgumentParser(description='A program to write a well-index map.')
  parser.add_argument('-i', '--index_max', required=True, help='Maximum index value.')
  parser.add_argument('-a', '--across_row_first', required=True, help='Index increases across row first (0=False; 1=True).')
  args = parser.parse_args()

  imax = int(args.index_max)
  across_row_first = bool(int(args.across_row_first))

  print('across_row_first: %s' % (across_row_first))

  print('  0  none     AAAAAAA')
  for i in range(imax):
    (ipl, well_id) = index_to_well(i, across_row_first)
    index_base_4 = base_convert(i+1, 4)
    encoded_index = index_base_4_to_encoded_index(index_base_4)
    print('%3d  P%02d-%s  %s' % (i+1, ipl+1, well_id, encoded_index))
