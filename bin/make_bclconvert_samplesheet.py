#!/usr/bin/env python3

# read:
#   samplesheet json file
#   p7 PCR sequence file
#   p5 PCR sequence file
#
# sort p7 PCR sequence file contents by well row and make a well to index map
# sort p5 PCR sequence file contents by well column and make a well to index map
#
# make bclconvert samplesheet with format
# where ip7 and ip5 are p7 and p5 well indices, respectively
#
# [Header]
# FileFormatVersion,2
# 
# [BCLConvert_Settings]
# CreateFastqForIndexReads,0
# BarcodeMismatchesIndex1, 1
# BarcodeMismatchesIndex2, 1
# 
# [BCLConvert_Data]
# Sample_ID,index,index2
# ip7_ip5,TTGGCAAGCC,TTGGTCATAA
# ip7_ip5,TTGGCAAGCC,TTGAGTATCA
# ip7_ip5,TTGGCAAGCC,GGTAGGTCTC
# 
# 
# Read rt_file, ligation_file, p7_file, and p5_file from json file.
# If any are zero length strings, using the barcode files in the
# pipeline repository. The paths are coded in this file.
#
# Functionality:
#   o  p7 and p5 barcodes
#      o  read files
#      o  convert well-sequence to index-sequence tables: use well-index conversion with input well-sequence table
#      o  make a dictionary keyed by lane with p7/p5 barcode paths as values.
#      o  allow for zero-length barcode file paths, in which case we use the (default) paths given as command line
#         parameters to this program.
#   o  samplesheet json file
#      o  read samplesheet file
#      o  extract p7 and p5 indices by lane
#      o  note: use the barcode indices for the sample_id in the bcl-convert samplesheet
#   o  make well-index dictionaries for the p7 and p5 barcodes
#   o  make list of output file entry data. Each entry consists of
#        o  lane
#        o  p7 seq (zero length string denotes dark)
#        o  p5 seq (zero length string denotes dark)
#        o  'sample_id' consists of p7 and p5 indices concatenated (0 index denotes dark)
#        o  either p7 or p5 may be 'dark' but not both

import sys
import json
import argparse
import re


#
# Program version string.
#
program_version = '0.2.0'


#
# ** Maximum number of microtiter plates.
#
MAX_NUM_PLATES = 8


#
# Read barcode file and store a dictionary keyed by well name.
# Standardize well names to upper case alphabetic characters and
# no plate name (for now).
#
def read_barcodes(filename):
  well_barcode_dict = {}
  try:
    with open(filename, 'r') as fh:
      for line in fh:
        toks = line.split()
        ntoks = len(toks)
        if(ntoks != 2):
          print('Error: read_barcodes: unexpected number of tokens: %d' % (ntoks))
          sys.exit(-1)
        well = toks[0].upper()
        seq  = toks[1]
        well_barcode_dict[well] = seq
  except:
    print('Error: unable to open barcode file \"%s\"' % (filename), sys.stderr)
    sys.exit(-1)
  return(well_barcode_dict)


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


#
# Note: The indexes begin with 1 so well A01 has index 1.
#
def make_well_to_index_dict(across_row_first):
  max_index = MAX_NUM_PLATES * 96
  well_index_dict = {}
  for index in range(1, max_index+1):
    (ipl, well_id) = index_to_well(index-1, across_row_first)
    well_str = 'P%02d-%s' % (ipl + 1, well_id)
    well_index_dict[well_str] = index
  return(well_index_dict)

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


#
# Get PCR data from samplesheet JSON file contents.
#
# json_data['sample_index_list']['lanes']
# json_data['sample_index_list']['ranges']
# json_data['sample_index_list']['p7_file']
# json_data['sample_index_list']['p5_file']
# 
def get_pcr_data(json_data):
  pcr_data_list= []
  for sample_index_dict in json_data['sample_index_list']:
    index_ranges = sample_index_dict['ranges'].split(':')
    p7_index_list = expand_index_list(index_ranges[1])
    p5_index_list = expand_index_list(index_ranges[2])
    lane_index_list = expand_index_list(sample_index_dict['lanes'])
    p7_file = sample_index_dict['p7_file']
    p5_file = sample_index_dict['p5_file']
    pcr_data_list.append( {'lane_index_list' : lane_index_list, 'p7_index_list' : p7_index_list, 'p5_index_list' : p5_index_list, 'p7_file' : p7_file, 'p5_file' : p5_file})
  return(pcr_data_list)


#
# Make a lane-to-p7+p5 barcode file path map using a dictionary
# keyed by the lane number with dictionary values with the p7 and p5
# barcode file paths.
#
def make_barcode_path_dict(pcr_data_list, default_p7_file, default_p5_file):
  barcode_path_dict = {}
  for sample_dict in pcr_data_list:
    p7_file = sample_dict['p7_file'] if len(sample_dict['p7_file']) > 0 else default_p7_file
    p5_file = sample_dict['p5_file'] if len(sample_dict['p5_file']) > 0 else default_p5_file
    for lane_index in sample_dict['lane_index_list']:
      lane_pcr_paths = barcode_path_dict.setdefault(lane_index, {'p7_file' : p7_file, 'p5_file' : p5_file })
      if(lane_pcr_paths['p7_file'] != p7_file or lane_pcr_paths['p5_file'] != p5_file):
        print('make_barcode_path_dict: error: inconsistent PCR file path(s) for lane %d' % (lane_index), file=sys.stderr)
        sys.exit(-1)
  return( barcode_path_dict)


#
# Reverse complement sequence 'seq'.
#
revcomp = str.maketrans('ACGTacgtRYMKrymkVBHDvbhd', 'TGCAtgcaYRKMyrkmBVDHbvdh')
def reverse_complement(seq):
    return seq.translate(revcomp)[::-1]


#
# Make a lane-to-barcode sequence map using a dictionary keyed by the
# the lane number with dictionary values that are dictionaries keyed
# by 'p7' and 'p5'. Those have values that are dictionaries keyed by
# index and the values are the corresponding sequences. So barcode_seq_dict
# is a dictionary of dictionaries of dictionaries. These look like
#
#   barcode_seq_dict[lane][(p7|p5)][index] = sequence
#
def make_barcode_seq_dict(barcode_path_dict, p7_well_to_index_dict, p5_well_to_index_dict, p5_rcmp):
  barcode_seq_dict = {}
  for lane in barcode_path_dict.keys():
    barcode_seq_dict[lane] = {}
    p7_file = barcode_path_dict[lane]['p7_file']
    p7_barcodes = read_barcodes(p7_file)
    barcode_seq_dict[lane]['p7'] = {}
    for well in p7_barcodes.keys():
      barcode_seq_dict[lane]['p7'][p7_well_to_index_dict[well]] = p7_barcodes[well]
    p5_file = barcode_path_dict[lane]['p5_file']
    p5_barcodes = read_barcodes(p5_file)
    barcode_seq_dict[lane]['p5'] = {}
    for well in p5_barcodes.keys():
      if(p5_rcmp == 'True'):
        barcode_seq_dict[lane]['p5'][p5_well_to_index_dict[well]] = reverse_complement(p5_barcodes[well])
      else:
        barcode_seq_dict[lane]['p5'][p5_well_to_index_dict[well]] = p5_barcodes[well]
  return(barcode_seq_dict)


def write_samplesheet_header(p7_index_zero, p5_index_zero, file_handle):
  print('[Header]', file=file_handle)
  print('FileFormatVersion,2', file=file_handle)
  print(file=file_handle)

  print('[BCLConvert_Settings]', file=file_handle)
  print('CreateFastqForIndexReads,0', file=file_handle)
  if(p7_index_zero == False):
    print('BarcodeMismatchesIndex1, 1', file=file_handle)
  if(p5_index_zero == False):
    print('BarcodeMismatchesIndex2, 1', file=file_handle)
  print(file=file_handle)
  return(0)


def write_samplesheet_data(distinct_tuples, barcode_seq_dict, p7_index_zero, p5_index_zero, file_handle):
  # Write data section header.
  print('[BCLConvert_Data]', file=file_handle)
  if(p7_index_zero == False and p5_index_zero == False):
    print('Lane,Sample_ID,index,index2', file=file_handle)
  elif(p7_index_zero == True):
    print('Lane,Sample_ID,index2', file=file_handle)
  elif(p5_index_zero == True):
    print('Lane,Sample_ID,index', file=file_handle)

  # Write data rows.
  for a_tuple in sorted(distinct_tuples):
    lane = a_tuple[0]
    p7_index = a_tuple[1]
    p5_index = a_tuple[2]
    if(p7_index != 0 and p5_index != 0):
      print('%d,%03d_%03d_%03d,%s,%s' % (lane,
                                         lane, p7_index, p5_index,
                                         barcode_seq_dict[lane]['p7'][p7_index],
                                         barcode_seq_dict[lane]['p5'][p5_index]), file=file_handle)
    elif(p7_index != 0):
      print('%d,%03d_%03d_%03d,%s' % (lane,
                                       lane, p7_index, p5_index,
                                       barcode_seq_dict[lane]['p7'][p7_index]), file=file_handle)
    elif(p5_index != 0):
      print('%d,%03d_%03d_%03d,%s' % (lane,
                                       lane, p7_index, p5_index,
                                       barcode_seq_dict[lane]['p5'][p5_index]), file=file_handle)


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


def check_distinct_tuples(distinct_tuples):
  # Check that all libraries specify p7 and p5 or p7 only or p5 only.
  p7_state = None
  p5_state = None
  for a_tuple in distinct_tuples:
    p7_index_zero = (a_tuple[1] == 0)
    p5_index_zero = (a_tuple[2] == 0)
    if(p7_state != None):
      if(p7_index_zero != p7_state):
        print('make_bclconvert_samplesheet: error: P7 primers include both none and non-none.', file=sys.stderr)
        sys.exit(-1)
    else:
      p7_state = p7_index_zero
    if(p5_state != None):
      if(p5_index_zero != p5_state):
        print('make_bclconvert_samplesheet: error: P7 primers include both none and non-none.', file=sys.stderr)
        sys.exit(-1)
    else:
      p5_state = p5_index_zero

  if(p7_index_zero == True and p5_index_zero == True):
    print('make_bclconvert_samplesheet: error: both p7 and p5 barcodes are set to zero.', file=sys.stderr)
    sys.exit(-1)
  return((p7_index_zero, p5_index_zero))


def  make_bclconvert_samplesheet(distinct_tuples, p7_index_zero, p5_index_zero, barcode_seq_dict, filename):
  with open(filename, 'w') as fh:
    write_samplesheet_header(p7_index_zero, p5_index_zero, fh)
    write_samplesheet_data(distinct_tuples, barcode_seq_dict, p7_index_zero, p5_index_zero, fh)


##############################
### Diagnostic code only.
##############################

# Consider filtering output using
#   grep '^diag' | awk '{print $3, $4, $5, $6, $7, $8}' | sort -k 1,2 -n | uniq
#

#
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


# 
#
#
def dump_pcr_well_barcodes(distinct_tuples, barcode_seq_dict):
  for tuple in distinct_tuples:
    lane = tuple[0]
    p7_index = tuple[1]
    p7_seq = barcode_seq_dict[lane]['p7'][p7_index]

    p5_index = tuple[2]
    if(p5_index > 0):
      p5_seq = barcode_seq_dict[lane]['p5'][p5_index]
    else:
      p5_seq = 'none'
    print('diag: %d %d %d %s %s %s %s' % (lane, p7_index, p5_index, index_to_well(p7_index-1, True)[1], p7_seq, index_to_well(p5_index-1, False)[1], p5_seq))

###
## End diagnostic code.
##############################


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to make bcl-convert samplesheet file.')
  parser.add_argument('-i', '--input', required=False, default=None, help='Input JSON samplesheet filename (required string).')
  parser.add_argument('-o', '--output', required=False, default=None, help='Output bcl-convert samplesheet filename (required string).')
  parser.add_argument('-7', '--p7_file', required=False, default=None, help='P7 barcode file path.')
  parser.add_argument('-5', '--p5_file', required=False, default=None, help='P5 barcode file path.')
  parser.add_argument('-r', '--p5_rcmp', required=True, default=None, help='P5 reverse complement (bool).')
  parser.add_argument('-v', '--version', required=False, default=None, help='Write version string to stdout.')
  args = parser.parse_args()

  # Write versions.
  if( args.version ):
    print( 'Program version: %s' % ( program_version ) )
    sys.exit( 0 ) 

  #
  # Make PCR barcode well to index dictionaries.
  #
  p7_well_to_index_dict = make_well_to_index_dict(True)
  p5_well_to_index_dict = make_well_to_index_dict(False)
#  print('p7 well to index dict: ', p7_well_to_index_dict)

  #
  # Read input samplesheet file.
  #
  json_data = read_json(args.input)
#  print(json.dumps(json_data['sample_index_list'], indent=2))
  pcr_data_list = get_pcr_data(json_data)

  #
  # Make lane-to-PCR barcode file path dictionary.
  #
  barcode_path_dict = make_barcode_path_dict(pcr_data_list, args.p7_file, args.p5_file)

  #
  # Make lane-to-barcode sequence dictionary.
  #
  barcode_seq_dict = make_barcode_seq_dict(barcode_path_dict, p7_well_to_index_dict, p5_well_to_index_dict, args.p5_rcmp)
#  print(json.dumps(barcode_seq_dict, indent=2))

  #
  # Make distinct index tuples of (<lane>, <p7_index>, <p5_index>).
  #
  distinct_tuples = make_distinct_tuples(pcr_data_list)

  #
  # Check tuples.
  #
  (p7_index_zero, p5_index_zero) = check_distinct_tuples(distinct_tuples)

  #
  # Make bcl-convert samplesheet file.
  #
  make_bclconvert_samplesheet(distinct_tuples, p7_index_zero, p5_index_zero, barcode_seq_dict, args.output)

  dump_pcr_well_barcodes(distinct_tuples, barcode_seq_dict)

