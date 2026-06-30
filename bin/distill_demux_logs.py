#!/usr/bin/env python3

import sys
import json
import argparse
import re

#
# Program version string.
#
program_version = '0.1.0'


#
# Maximum number of lanes/wafers.
#
max_lanes = 16


#
# Maximum well index.
#
max_plate = 16
max_well_index = max_plate * 96


#
# Read log json file.
#
def read_json(filename):
  with open(filename, 'r') as fh:
    json_data = json.load(fh)
  return(json_data)


def initialize_pcr_counts_dict():
  pcr_counts_dict = dict()
  pcr_counts_dict['input_read_pairs'] = 0
  pcr_counts_dict['rt_9_matched'] = 0
  pcr_counts_dict['lg_9_matched'] = 0
  pcr_counts_dict['rtlg_9_matched'] = 0
  pcr_counts_dict['rt_10_matched'] = 0
  pcr_counts_dict['lg_10_matched'] = 0
  pcr_counts_dict['rtlg_10_matched'] = 0
  pcr_counts_dict['rtlg_9_and_rtlg_10_matched'] = 0
  pcr_counts_dict['matches_failed'] = 0
  return(pcr_counts_dict)


def sum_pcr_counts(json_data, pcr_counts_dict):
  pcr_counts_dict['input_read_pairs'] += json_data['input_read_pairs']
  pcr_counts_dict['rt_9_matched'] += json_data['rt_9_matched']
  pcr_counts_dict['lg_9_matched'] += json_data['lg_9_matched']
  pcr_counts_dict['rtlg_9_matched'] += json_data['rtlg_9_matched']
  pcr_counts_dict['rt_10_matched'] += json_data['rt_10_matched']
  pcr_counts_dict['lg_10_matched'] += json_data['lg_10_matched']
  pcr_counts_dict['rtlg_10_matched'] += json_data['rtlg_10_matched']
  pcr_counts_dict['rtlg_9_and_rtlg_10_matched'] += json_data['rtlg_9_and_rtlg_10_matched']
  pcr_counts_dict['matches_failed'] += json_data['matches_failed']
  return(0)


def write_pcr_counts(pcr_counts_dict, ofp):
  print('input_read_pairs:           %d' % (pcr_counts_dict['input_read_pairs']), file=ofp)
  print('rt_9_matched:               %d (%.03f)' % (pcr_counts_dict['rt_9_matched'], pcr_counts_dict['rt_9_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('lg_9_matched:               %d (%.03f)' % (pcr_counts_dict['lg_9_matched'], pcr_counts_dict['lg_9_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('rtlg_9_matched:             %d (%.03f)' % (pcr_counts_dict['rtlg_9_matched'], pcr_counts_dict['rtlg_9_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('rt_10_matched:              %d (%.03f)' % (pcr_counts_dict['rt_10_matched'], pcr_counts_dict['rt_10_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('lg_10_matched:              %d (%.03f)' % (pcr_counts_dict['lg_10_matched'], pcr_counts_dict['lg_10_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('rtlg_10_matched:            %d (%.03f)' % (pcr_counts_dict['rtlg_10_matched'], pcr_counts_dict['rtlg_10_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('rtlg_9_and_rtlg_10_matched: %d (%.03f)' % (pcr_counts_dict['rtlg_9_and_rtlg_10_matched'], pcr_counts_dict['rtlg_9_and_rtlg_10_matched'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print('matches_succeeded:          %d (%.03f)' % (pcr_counts_dict['input_read_pairs'] - pcr_counts_dict['matches_failed'],
                                                    1 - (pcr_counts_dict['matches_failed'] / pcr_counts_dict['input_read_pairs'])), file=ofp)
  print('matches_failed:             %d (%.03f)' % (pcr_counts_dict['matches_failed'], pcr_counts_dict['matches_failed'] / pcr_counts_dict['input_read_pairs']), file=ofp)
  print()

  return(0)


def initialize_well_counts_dict():
  rt_count_list = [0] * max_well_index
  well_counts_dict = {'rt_counts': rt_count_list}
  return(well_counts_dict)


def sum_well_counts(json_data, well_counts_dict):
  rt_counter_vec = json_data['rt_counter_vec']
  for rt_count_dict in rt_counter_vec:
    counts = rt_count_dict['whitelist_read_counts'] + rt_count_dict['mismatch_read_counts']
    rt_well_index = rt_count_dict['well_index']
    well_counts_dict['rt_counts'][rt_well_index] += counts
  return(0)


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


# OK
# Note: the first true well_index value is '0'. The
#       well_index value < 0 represents no well.
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


def write_well_counts(well_counts_dict, ofp):
  rt_count_list = well_counts_dict['rt_counts']
  total_counts = 0
  for rt_well_count in rt_count_list:
    total_counts += rt_well_count  
  print('rt well index counts (%d total)' % (total_counts))
  for idx, rt_well_count in enumerate(rt_count_list):
    if(rt_well_count > 0 ):
      well_tuple = index_to_well(idx-1, True)
      swell = 'P%02d-%s' % (well_tuple[0]+1, well_tuple[1])
      print('%4d  %s  %d (%.03f)' % (idx, swell, rt_well_count, rt_well_count / total_counts), file=ofp)
  print()
  return(0)


def initialize_marginal_counts_dict():
  lane_count_list = [0] * max_lanes
  p7_count_list   = [0] * (max_well_index + 1)
  p5_count_list   = [0] * (max_well_index + 1)
  marginal_counts_dict = {'lane_count_list': lane_count_list, 'p7_count_list': p7_count_list, 'p5_count_list': p5_count_list}
  return(marginal_counts_dict)


def sum_marginal_counts(json_data, marginal_counts_dict):
  lane_count_list = marginal_counts_dict['lane_count_list']
  p7_count_list   = marginal_counts_dict['p7_count_list']
  p5_count_list   = marginal_counts_dict['p5_count_list']

  rt_counter_vec = json_data['rt_counter_vec']
  for rt_count_dict in rt_counter_vec:
    lane_index = rt_count_dict['lanes'][0]
    p7_index   = rt_count_dict['p7_indices'][0]
    p5_index   = rt_count_dict['p5_indices'][0]
    count      = rt_count_dict['whitelist_read_counts'] + rt_count_dict['mismatch_read_counts']

    lane_count_list[lane_index] += count
    p7_count_list[p7_index] += count
    p5_count_list[p5_index] += count

  return(0)


def write_marginal_counts(marginal_counts_dict, ofp):
  lane_count_list = marginal_counts_dict['lane_count_list']
  p7_count_list   = marginal_counts_dict['p7_count_list']
  p5_count_list   = marginal_counts_dict['p5_count_list']

  total_counts = 0
  for lane_count in lane_count_list:
    total_counts += lane_count
  print('lane marginal read counts (%d total)' % (total_counts))
  for idx, lane_count in enumerate(lane_count_list):
    if(lane_count > 0):
      print('%4d  %d  (%.03f)' % (idx, lane_count, lane_count / total_counts), file=ofp)
  print()

  total_counts = 0
  for p7_count in p7_count_list:
    total_counts += p7_count
  print('p7 marginal read counts (%d total)' % (total_counts))
  for idx, p7_count in enumerate(p7_count_list):
    if(p7_count > 0):
      print('%4d  %d  (%.03f)' % (idx, p7_count, p7_count / total_counts), file=ofp)
  print()

  total_counts = 0
  for p5_count in p5_count_list:
    total_counts += p5_count
  print('p5 marginal read counts (%d total)' % (total_counts))
  for idx, p5_count in enumerate(p5_count_list):
    if(p5_count > 0):
      print('%4d  %d  (%.03f)' % (idx, p5_count, p5_count / total_counts), file=ofp)
  print()


def write_all_counts_json(pcr_counts_dict, well_counts_dict, marginal_counts_dict):
  #
  # Store the counts in a dictionary.
  #
  counts_dict = dict()
  counts_dict['pcr_counts'] = pcr_counts_dict

  well_count_json = list()
  rt_count_list = well_counts_dict['rt_counts']
  for idx , rt_well_count in enumerate(rt_count_list):
    if(rt_well_count > 0 ):
      well_tuple = index_to_well(idx-1, True)
      swell = 'P%02d-%s' % (well_tuple[0]+1, well_tuple[1])
      count_list = [idx, swell, rt_well_count]
      well_count_json.append(count_list)
  counts_dict['rt_well_count_list'] = well_count_json

  lane_count_list = marginal_counts_dict['lane_count_list']
  lane_count_json = list()
  for idx , lane_count in enumerate(lane_count_list):
    if(lane_count > 0 ):
      count_list = [idx, lane_count]
      lane_count_json.append(count_list)
  counts_dict['lane_marginal_count_list'] = lane_count_json

  p7_count_list = marginal_counts_dict['p7_count_list']
  p7_count_json = list()
  for idx , p7_count in enumerate(p7_count_list):
    if(p7_count > 0 ):
      well_tuple = index_to_well(idx-1, True)
      swell = 'P%02d-%s' % (well_tuple[0]+1, well_tuple[1])
      count_list = [idx, swell, p7_count]
      p7_count_json.append(count_list)
  counts_dict['p7_well_marginal_count_list'] = p7_count_json
 
  p5_count_list = marginal_counts_dict['p5_count_list']
  p5_count_json = list()
  for idx , p5_count in enumerate(p5_count_list):
    if(p5_count > 0 ):
      well_tuple = index_to_well(idx-1, False)
      swell = 'P%02d-%s' % (well_tuple[0]+1, well_tuple[1])
      count_list = [idx, swell, p5_count]
      p5_count_json.append(count_list)
  counts_dict['p5_well_marginal_count_list'] = p5_count_json

  print('counts dict:')
  print(json.dumps(counts_dict, indent=2))

  return(0)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to distill demux JSON log files.')
  parser.add_argument('-i', '--input', required=True, default=None, nargs='+', help='Input JSON log filenames (required string(s)).')
#  parser.add_argument('-o', '--output', required=False, default=None, help='Output JSON filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Initialize counters.
  #
  pcr_counts_dict  = initialize_pcr_counts_dict()
  well_counts_dict = initialize_well_counts_dict()
  marginal_counts_dict = initialize_marginal_counts_dict()

  #
  # Read log JSON files and make counts.
  #
  for filename in args.input:
    json_data = read_json(filename)
    sum_pcr_counts(json_data, pcr_counts_dict)
    sum_well_counts(json_data, well_counts_dict)
    sum_marginal_counts(json_data, marginal_counts_dict)

  #
  # Write counts.
  #
#  write_pcr_counts(pcr_counts_dict, sys.stdout)
#  write_well_counts(well_counts_dict, sys.stdout)
#  write_marginal_counts(marginal_counts_dict, sys.stdout)
  write_all_counts_json(pcr_counts_dict, well_counts_dict, marginal_counts_dict)

