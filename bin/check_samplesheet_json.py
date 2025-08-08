#!/usr/bin/env python3

import sys
import json
import argparse
import re


program_version = '0.1.0'

#
# Check JSON samplesheet file.
#   o  check for duplicate combinations of process_group, lanes, rt indices, p7 indices, and p5 indices (done)
#   o  check for duplicate sample names in each process group (I don't see how to do this.)
#   o  check for different hash files in a process group
#   o  check for different rt files in a process group
#   o  check for different ligation files in a process group
#   o  check for different p7 files in a process group
#   o  check for different p5 files in a process group
#   o  check for different p7 files in a lane
#   o  check for different p5 files in a lane
#


#  sample_index_list.append( { 'sample_id' : row_out['sample_name'],
#                              'ranges' : ':'.join( [ make_index_string( row_out['rt_index_list'] ),
#                                                     make_index_string( row_out['p7_index_list'] ),
#                                                     make_index_string( row_out['p5_index_list'] ) ] ),
#                              'lanes' : lanes,
#                              'tissue' : row_out['tissue'],
#                              'genome' : row_out['genome'],
#                              'hash_file' : row_out['hash_file'],
#                              'sample_flags' : ','.join( row_out['sample_flags'] ),
#                              'external_sample_name' : row_out['external_sample_name'],
#                              'wrap_group' : row_out['wrap_group'],
#                              'rt_file' : row_out['rt_file'],
#                              'ligation_file' : row_out['ligation_file'],
#                              'p7_file' : row_out['p7_file'],
#                              'p5_file' : row_out['p5_file'],
#                              'library' : row_out['library'],
#                              'process_group' : process_group })
#

#
# Read samplesheet json file.
#
def read_json(filename):
  with open(filename, 'r') as fh:
    json_data = json.load(fh)
  return(json_data)


#
# Expand index list.
# Example: 1-4,8-12
# Return a list of distinct indices.
regex_pattern = r'([0-9]+)([-]([0-9]+))?$'
def expand_index_list(index_string):
  index_list = []
  for index_spec in index_string.split(','):
    mobj = re.match(regex_pattern, index_spec)
    if(mobj == None):
      print('Error: expand_index_list: bad index specification: %s' % (index_spec))
      sys.exit(-1)
    index1 = int(mobj.group(1))
    index2 = index1
    if(mobj.group(2) != None):
      index2 = int(mobj.group(3))
    for one_index in range(index1, index2+1):
      index_list.append(one_index)
  return(list(set(index_list)))



def make_distinct_tuples(pcr_data_list):
  # Make distinct (lane, p7_index, p5_index) tuples.
  sample_tuples = []
  for sample_dict in pcr_data_list:
    lane_index_list = sample_dict['lane_index_list']
    p7_index_list = sample_dict['p7_index_list']
    p5_index_list = sample_dict['p5_index_list']
    for lane in lane_index_list:
      for p7_index in p7_index_list:
        for p5_index in p5_index_list:
          sample_tuples.append((lane, p7_index, p5_index))
  distinct_tuples = set(sample_tuples)
  return(distinct_tuples)


#
# Make an index list from a string of index ranges.
#
index_range_re_pattern = r'([0-9]+)([-]([0-9]+))?$'
def make_index_list( index_string ):
  index_list = []
  for index_range in index_string.split( ',' ):
    mobj = re.match( index_range_re_pattern, index_range.strip() )
    index1 = int( mobj.group( 1 ) )
    index2 = index1
    if( mobj.group( 2 ) ):
      index2 = int( mobj.group( 3 ) )
    for i in range( index1, index2 + 1 ):
      index_list.append( i )
  return( index_list )


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


def pad_well_col(well_col, zero_pad, id_length):
  if zero_pad:
      template = '%%0%sd' % id_length
  else:
      template = '%s'
  col_id = template % (well_col)
  return col_id


def index_string_to_well_string( index_string, across_row_first, show_plate ):
  well_string = ''
  for item in index_string.split(','):
    mobj = re.match( r'^([0-9]+)([-]([0-9]+))?$', item.strip() )
    if( mobj == None):
      print('Error: index_string_to_well_string: unable to parse index string', file=sys.stderr)
      sys.exit(-1)
    ipl, well = index_to_well(int( mobj.group( 1 ) ) - 1, across_row_first)
    if( ipl == 0 and not show_plate ):
      if( len( well_string ) > 0 ):
        well_string += ','
      well_string += well
    else:
      if( len( well_string ) > 0 ):
        well_string += ','
      well_string += 'P%02d-%s' % ( ipl + 1, well )
    if( mobj.group( 3 ) != None ):
      ipl, well = index_to_well(int( mobj.group( 3 ) ) - 1, across_row_first)
      if( ipl == 0 and not show_plate ):
        well_string += ':%s' % ( well )
      else:
        well_string += ':P%02d-%s' % ( ipl + 1, well )

  return(well_string)


def make_index_string( index_list ):
  """
  Convert a list of (integer) barcode well indices to an index string where
    o  repeated indices are dropped; that is, keep only distinct indices
    o  sequences of counting numbers are expressed as ranges, for example, 5 6 7 8 9 => 5-9
    o  indices and index ranges are separated by commas
  """
  index_string = ''
  if( len( index_list ) == 0 ):
    return( index_string )
  index_list.sort()
  index_prev = None
  index1 = None
  for i in index_list:
    if( index_prev ):
      if( i == index_prev ):
        continue
      elif( i > index_prev + 1 ):
        if( len( index_string ) > 0 ):
          index_string += ','
        if( index_prev > index1 ):
          index_string += '%d-%d' % ( index1, index_prev )
        else:
          index_string += '%d' % ( index_prev )
        index1 = i
    else:
      index1 = i
    index_prev = i
  # last index in list
  if( len( index_string ) > 0 ):
    index_string += ','
  if( index_prev > index1 ):
    index_string += '%d-%d' % ( index1, index_prev )
  else:
    index_string += '%d' % ( index_prev )
  return( index_string )


def write_sample_index_dict(sample_index_dict, prefix):
  print( '%ssample name: %s' % (prefix, sample_index_dict['sample_id']))
  print( '%sgenome: %s' % (prefix, sample_index_dict['genome']))
  print( '%sprocess group: %s' % (prefix, sample_index_dict['process_group']))
  print( '%slanes: %s' % (prefix, sample_index_dict['lanes']))
  (rt_indices, p7_indices, p5_indices) = sample_index_dict['ranges'].split(':')
  print( '%srt wells: %s' % (prefix, index_string_to_well_string(rt_indices, True, True)))
  print( '%sp7 wells: %s' % (prefix, index_string_to_well_string(p7_indices, True, True)))
  print( '%sp5 wells: %s' % (prefix, index_string_to_well_string(p5_indices, False, True)))
  if(len(sample_index_dict['hash_file']) != 0):
    print( '%shash file: %s' % (prefix, sample_index_dict['hash_file']))
  if(len(sample_index_dict['rt_file']) != 0):
    print( '%srt file: %s' % (prefix, sample_index_dict['rt_file']))
  if(len(sample_index_dict['p7_file']) != 0):
    print( '%sp7 file: %s' % (prefix, sample_index_dict['p7_file']))
  if(len(sample_index_dict['p5_file']) != 0):
    print( '%sp5 file: %s' % (prefix, sample_index_dict['p5_file']))
  if(len(sample_index_dict['ligation_file']) != 0):
    print( '%sligation file: %s' % (prefix, sample_index_dict['ligation_file']))

  return( 0 )


#
# Check that no rows have the same set of
#   o  process group
#   o  lanes
#   o  rt indices
#   o  p7 indices
#   o  p5 indices
#
# Strategy
#   o  form strings that store the information
#   o  expand the lanes, rt, ligation, p7, and p5 indices to check for intersections
#   o  sort the strings
#   o  check for duplicate strings
#
def check_for_duplicate_combinations(sample_index_list):
  error_flag = False
  # Use sample_index_subdict_lookup for reporting duplicates.
  sample_index_subdict_lookup = dict()
  entry_list = list()
  for sample_index_dict in sample_index_list:
    process_group = sample_index_dict['process_group']
    (rt_indices, p7_indices, p5_indices) = sample_index_dict['ranges'].split(':')
    sample_index_subdict = {'sample_id' : sample_index_dict['sample_id'],
                            'genome' : sample_index_dict['genome'],
                            'process_group' : process_group,
                            'lanes' : sample_index_dict['lanes'],
                            'ranges' : sample_index_dict['ranges'],
                            'hash_file' : sample_index_dict['hash_file'],
                            'rt_file' : sample_index_dict['rt_file'],
                            'ligation_file' : sample_index_dict['ligation_file'],
                            'p7_file' : sample_index_dict['p7_file'],
                            'p5_file' : sample_index_dict['p5_file']}
    lanes_list = make_index_list(sample_index_dict['lanes'])
    rt_index_list = make_index_list(rt_indices)
    p7_index_list = make_index_list(p7_indices)
    p5_index_list = make_index_list(p5_indices)
    for lane in lanes_list:
      for rt_index in rt_index_list:
        for p7_index in p7_index_list:
          for p5_index in p5_index_list:
            entry_string = '%s|%s|%s|%s|%s' % (process_group, 
                                               lane,
                                               rt_index,
                                               p7_index,
                                               p5_index)
            entry_list.append(entry_string)
            if(not entry_string in sample_index_subdict_lookup):
              sample_index_subdict_lookup[entry_string] = list()
            sample_index_subdict_lookup[entry_string].append(sample_index_subdict)
  entry_list.sort()
  dup_entry_list = list()
  for i in range(1,len(entry_list)):
    if(entry_list[i] == entry_list[i-1]):
      if(not entry_list[i] in dup_entry_list):
        dup_entry_list.append(entry_list[i])
  if(len(dup_entry_list) != 0):
    error_flag = True
    dup_subdict_list = []
    for dup_entry in dup_entry_list:
      for subdict in sample_index_subdict_lookup[dup_entry]:
        if(not subdict in dup_subdict_list):
          dup_subdict_list.append(subdict)
    print('Error: duplicate combination of process_group, lanes, rt_indices, p7_indices, and p5_indices')
    for subdict in dup_subdict_list:
      write_sample_index_dict(subdict, '  ')
      print('  --')
  return(error_flag)


def check_file_process_group(sample_index_list, file_type):
  error_flag = False
  file_dict = dict()
  for sample_index_dict in sample_index_list:
    process_group = sample_index_dict['process_group']
    file_name = sample_index_dict[file_type]  
    if(not process_group in file_dict):
      file_dict[process_group] = list()
    if(not file_name in file_dict[process_group]):
      file_dict[process_group].append(file_name)

  for process_group in file_dict:
    if(len(file_dict[process_group]) > 1):
      error_flag = True
      print('Error: more than one %s in process group \'%s\'' % (file_type, process_group))
      for file_name in file_dict[process_group]:
        if(len(file_name) > 0):
          print('  %s' % (file_name))
        else:
          print('  <empty string>')
  return(error_flag)


def check_file_lanes(sample_index_list, file_type):
  error_flag = False
  file_dict = dict()
  for sample_index_dict in sample_index_list:
    lanes = sample_index_dict['lanes']
    lanes_list = make_index_list(lanes)
    file_name = sample_index_dict[file_type]
    for lane in lanes_list:
      if(not lane in file_dict):
        file_dict[lane] = list()
      if(not file_name in file_dict[lane]):
        file_dict[lane].append(file_name)

  for lane in file_dict:
    if(len(file_dict[lane]) > 1):
      error_flag = True
      print('Error: more than one %s in lane \'%s\'' % (file_type, lane))
      for file_name in file_dict[lane]:
        if(len(file_name) > 0):
          print('  %s' % (file_name))
        else:
          print('  <empty string>')
  return(error_flag)


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to JSON samplesheet file.')
  parser.add_argument('-i', '--input', required=True, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  
  json_data = read_json(args.input)
  sample_index_list = json_data['sample_index_list']

  exit_flag = 0
  error_flag = check_for_duplicate_combinations(sample_index_list)
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_process_group(sample_index_list, 'hash_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_process_group(sample_index_list, 'rt_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_lanes(sample_index_list, 'rt_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_process_group(sample_index_list, 'ligation_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_lanes(sample_index_list, 'ligation_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_process_group(sample_index_list, 'p7_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_lanes(sample_index_list, 'p7_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_process_group(sample_index_list, 'p5_file')
  if(error_flag):
    exit_flag = 1

  error_flag = check_file_lanes(sample_index_list, 'p5_file')
  if(error_flag):
    exit_flag = 1

  sys.exit(exit_flag)

