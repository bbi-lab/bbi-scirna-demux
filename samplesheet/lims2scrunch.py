#!/usr/bin/env python3

#
# Convert a LIMS CSV manifest file to a bbi-dmux samplesheet file.
# The manifest file format is set on 20231025. It may change in
# the future.
#
# This a minimal script at this time.
#
# Notes:
#   o  the values in the parameter.txt file are applied to all
#      samples in the output file, except for the *_sample_name
#      values, which are applied only to the specified samples.
#      Consequently, you may need to edit some of the values in
#      the output .csv file using a spreadsheet program.
#


import sys
import argparse
import csv
import json
import re
import operator
import shlex


#
# This is a multi-line description string for argparse help.
#
description_string = '''Convert a BBI LIMS manifest to a sample-oriented CSV samplesheet file.

Add columns to the output file using a parameters input file, which has values common to all rows. The template for the parameters file is

# rt_file: ""
# ligation_file: ""
# p7_file: ""
# p5_file: ""
# hash_file: ""
# lanes: ""
# p7_wells: ""
# p5_wells: ""
# barnyard_sample_name: ""
# sentinel_sample_name: ""
# keyhole_sample_name: ""

To set a parameter, remove the '#' and add the value between the pair of double quotes. This script uses the shlex python package to interpret the parameter file contents.

You will almost certainly need to edit the resulting CSV file using a spreadsheet program before running the scirna_samplesheet.py script on it in order to make the samplesheet JSON file required by the bbi-scirna-demux and bbi-scirna-analyze pipelines.

Specify lanes individually, separated by commas, or as ranges, using a colon between the first and last lanes in the range. As examples, "1,3,6" or "1:8" or "2:4,6:8".

Specify p7_wells and p5_wells individually, separated by commas, or as ranges, using a colon between the first and last wells.

Note that the p7 well ranges increase along rows so "P01-A01:P01-A12" is the first p7 row, and the p5 well ranges increase along columns so "P01-A01:P01-H01" is the first p5 column. An entire plate is "P01-A01:P01-H12".

These pipelines require that the p7 and p5 well names include the plate number.
'''



# "name","experiment","assay","plate","coordinates","RT_block","BBI_ID","investigatorSpecimenId","investigator","organism","tissue","genome"
# "RNA3-074-P01-A1-24.0262","RNA3-074","sci-RNA-seq","P01","A1","RNA3-074_24.0262","24.0262","RA-gastruloid-2day","Hamazaki","Human","Gastruloid","Human"
# "RNA3-074-P01-A10-24.0270","RNA3-074","sci-RNA-seq","P01","A10","RNA3-074_24.0270","24.0270","RA-gastruloid-36h","Hamazaki","Human","Gastruloid","Human"
# "RNA3-074-P01-A11-24.0270","RNA3-074","sci-RNA-seq","P01","A11","RNA3-074_24.0270","24.0270","RA-gastruloid-36h","Hamazaki","Human","Gastruloid","Human"

# Check the LIMS manifest CSV file header.
def check_header(row):
  if(row[0] !=  'name' or
     row[1] !=  'experiment' or
     row[2] !=  'assay' or
     row[3] !=  'plate' or
     row[4] !=  'coordinates' or
     row[5] !=  'RT_block' or
     row[6] !=  'BBI_ID' or
     row[7] !=  'investigatorSpecimenId' or
     row[8] !=  'investigator' or
     row[9] !=  'organism' or
     row[10] != 'tissue' or
     row[11] != 'genome'):
    print('Error: unexpected header token in row: ', row, file=sys.stderr)
    sys.exit(1)


# Read the LIMS manifest CSV file.
def read_lims_manifest(infile):
  rows = []
  with open(infile, newline='') as ifp:
    reader = csv.reader(ifp, dialect='unix')
    # Read rows from input file.
    for irow, row in enumerate(reader):
      if(irow == 0):
        check_header(row)
        header_row = row
        continue
      rows.append(row)
  return((header_row, rows))


def store_row_values_dicts(header_row, inrows):
  pobj = re.compile(r'([A-H])([0-9]+)')
  inrows_dicts = []
  for inrow in inrows:
    dict_tmp = dict()
    for icol, cell_value in enumerate(inrow):
      dict_tmp[header_row[icol]] = cell_value
    plate = dict_tmp['plate']
    well_coordinate = dict_tmp['coordinates']
    mobj = pobj.match(well_coordinate)
    if(mobj == None):
      print('Unable to match well coordinate \'%s\'' % (well_coordinate), file=sys.stderr)
      sys.exit(1)
    dict_tmp['rt_well'] = '%s-%s%02d' % (plate, mobj.group(1), int(mobj.group(2)))
    inrows_dicts.append(dict_tmp)
  return(inrows_dicts)


def gather_sample_rows(inrows_dicts):
  sample_dicts = dict()
  for row_dict in inrows_dicts:
    rt_block = row_dict['RT_block']
    if(sample_dicts.get(rt_block) != None):
      if(row_dict['investigatorSpecimenId'] != sample_dicts[rt_block]['investigatorSpecimenId']):
        print('Error: inconsistent investigatorSpecimenId for sample %s' % (rt_block), file=sys.stderr)
      if(row_dict['investigator'] != sample_dicts[rt_block]['investigator']):
        print('Error: inconsistent investigatorSpecimenId for sample %s' % (rt_block), file=sys.stderr)
      if(row_dict['organism'] != sample_dicts[rt_block]['organism']):
        print('Error: inconsistent investigatorSpecimenId for sample %s' % (rt_block), file=sys.stderr)
      if(row_dict['tissue'] != sample_dicts[rt_block]['tissue']):
        print('Error: inconsistent investigatorSpecimenId for sample %s' % (rt_block), file=sys.stderr)
      if(row_dict['genome'] != sample_dicts[rt_block]['genome']):
        print('Error: inconsistent investigatorSpecimenId for sample %s' % (rt_block), file=sys.stderr)
      #
      # Append well.
      #
      sample_dicts[rt_block]['rt_wells'].append(row_dict['rt_well'])
    else:
      sample_dicts[rt_block] = dict()
      sample_dicts[rt_block]['investigatorSpecimenId'] = row_dict['investigatorSpecimenId']
      sample_dicts[rt_block]['investigator'] = row_dict['investigator']
      sample_dicts[rt_block]['organism'] = row_dict['organism']
      sample_dicts[rt_block]['tissue'] = row_dict['tissue']
      sample_dicts[rt_block]['genome'] = row_dict['genome']
      sample_dicts[rt_block]['rt_wells'] = list()
      sample_dicts[rt_block]['rt_wells'].append(row_dict['rt_well'])
  return(sample_dicts)
        


# Write samplesheet CSV file to outfile.
def write_sample_sheet(sample_dicts, parameter_dict, outfile):
  with open(outfile, 'w+', newline='') as ofp:
    writer = csv.writer(ofp, dialect='unix', quotechar='"', delimiter=',', escapechar='\\')
    writer.writerow(['sample_name',
                     'genome',
                     'process_group',
                     'lanes',
                     'rt_wells',
                     'p7_wells',
                     'p5_wells',
                     'rt_file',
                     'ligation_file',
                     'p7_file',
                     'p5_file',
                     'hash_file',
                     'sample_flags',
                     'tissue',
                     'external_sample_name',
                     'wrap_group'])
    for rt_block in sample_dicts.keys():
      rt_wells_string = ','.join(sample_dicts[rt_block]['rt_wells'])
      sample_flags = ''
      if(len(parameter_dict['barnyard_sample_name']) > 0 and rt_block == parameter_dict['barnyard_sample_name']):
        sample_flags = sample_flags + 'B'
      if(len(parameter_dict['sentinel_sample_name']) > 0 and rt_block == parameter_dict['sentinel_sample_name']):
        sample_flags = sample_flags + 'S'
      if(len(parameter_dict['keyhole_sample_name']) > 0 and rt_block == parameter_dict['keyhole_sample_name']):
        sample_flags = sample_flags + 'K'

      writer.writerow([rt_block,
                      sample_dicts[rt_block]['genome'],
                      '1',
                      parameter_dict['lanes'],
                      rt_wells_string,
                      parameter_dict['p7_wells'],
                      parameter_dict['p5_wells'],
                      parameter_dict['rt_file'],
                      parameter_dict['ligation_file'],
                      parameter_dict['p7_file'],
                      parameter_dict['p5_file'],
                      parameter_dict['hash_file'],
                      sample_flags,
                      sample_dicts[rt_block]['tissue'],
                      sample_dicts[rt_block]['investigatorSpecimenId'],
                      sample_dicts[rt_block]['investigator']])
  return(0)


# rt_file: ""
# lig_file: ""
# p7_file: ""
# p5_file: ""
# hash_file: ""
# lanes: ""
# p7_wells: ""
# p5_wells: ""
# lanes: ""
def read_parameter_file(filename):
  parameter_list = ['rt_file', 'ligation_file', 'p7_file', 'p5_file', 'hash_file', 'lanes', 'p7_wells', 'p5_wells', 'barnyard_sample_name', 'sentinel_sample_name', 'keyhole_sample_name']
  parameter_dict = dict()
  for parameter in parameter_list:
    parameter_dict[parameter] = ''
  print('filename: %s' % (filename))
  with open(filename, 'r') as fp:
    for line in fp:
      line_clean = line.strip()
      if(re.match('#', line_clean) or len(line_clean) == 0):
        continue
      parts = shlex.split(line)
      key = re.sub(r':', '', parts[0])
      if(key in parameter_list):
        parameter_dict[key] = parts[1]
      else:
        print('Error: unrecognized parameter (%s) in parameter file %s' % (key, filename), file=sys.stderr)
        sys.exit(1)
  print('parameter_dict: ')
  print(parameter_dict)
  return(parameter_dict)



if __name__ == '__main__':

  parser = argparse.ArgumentParser(description=description_string, formatter_class=argparse.RawDescriptionHelpFormatter)
  parser.add_argument('-i', '--in_file', required=True, help='Input LIMS manifest CSV filename.')
  parser.add_argument('-o', '--out_file', required=True, help='Output CSV filename.')
  parser.add_argument('-p', '--parameter_file', required=True, help='Input parameters filename.')

  args = parser.parse_args()

  infile = args.in_file
  outfile = args.out_file
  parameter_file = args.parameter_file

  print('Input file: %s' % (infile))
  print('Output file: %s' % (outfile))

  if(parameter_file != None):
    parameter_dict = read_parameter_file(parameter_file)
  else:
    parameter_dict = None

  # Read input file.
  (header_row, inrows) = read_lims_manifest(infile)

  # Store values in dictionaries.
  inrows_dicts = store_row_values_dicts(header_row, inrows)

  # Gather samples into single records.
  sample_dicts = gather_sample_rows(inrows_dicts)


  #
  # Write sample-oriented CSV file.
  #
  write_sample_sheet(sample_dicts, parameter_dict, outfile)


  # sample_dicts = gather_sample_rows(inrows_dicts)
  # print('%s' % (json.dumps(sample_dicts, indent=2)))


