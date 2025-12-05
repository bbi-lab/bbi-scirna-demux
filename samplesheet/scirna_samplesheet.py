#!/usr/bin/env python3

#
# TODO:
#   o  add process_group documentation
#
#   o  column_header_name_list = [ 'rt_wells', 'sample_name', 'genome', 'sample_flags', 'external_sample_name', 'tissue', 'wrap_group', 'lanes', 'hash_file', 'rt_file', 'ligation_file', 'p7_file', 'p5_file', 'library', 'process_group' ]
#   o  columns_required_list = [ 'rt', 'p5', 'p7', 'sample_name', 'genome' ]
#   o  check --lanes as optional parameter (uses -n <n> to define default)   (no)
#   o  point out that different samples may have different rt_file,
#      ligation_file, p7_file, p5_file, and hash_file values in case of
#      different libraries
#   o  check sample_flags description for completeness
#        o  the flags are used to identify certain sample types, which receive
#           type-specific processing.
#        o  the sample
#   o  explain how values are used in the pipelines
#


#
# Notes:
#   o  the master file is bbi-scirna-demux/samplesheet/scirna_samplesheet.py
#   o  use these bbi-scirna-* pipelines with caution. They have not been tested
#      thoroughly. For example, in cases where different lanes use different
#      rt, lig, p7, and/or p5 barcode files.
#

"""
Program: scirna_samplesheet.py
Summary:
  This program reads a (front-end) CSV spreadsheet file, performs a variety of
  checks, 'fixes' sample names, and writes a 'back-end' samplesheet file for
  use with the processing pipeline where the 'back-end' samplesheet file is
  in JSON format.

Input (front-end) samplesheet format:
  o  the input samplesheet file is a CSV format spreadsheet file (use a
     spreadsheet program to create the CSV file)
  o  the CSV column order is arbitrary (but must be consistent within the file)
  o  the header values do not depend on case
  o  checks rows:
       o  trims off extra cells from rows when a data row has more cells
          than the header line.
       o  ignores rows with all empty cells.
  o  the first row is a header with the following required columns:
       o  rt barcode identifier with value: 'rt_wells'.
       o  p7 barcode identifier with value: 'p7_wells' or 'p7_rows'.
       o  p5 barcode identifier with value: 'p5_wells' or 'p7_columns'.
       o  sample name identifier with value: 'sample_name'.
       o  genome label identifier with value: 'genome'.
     and the following optional column header values and column contents:
       o  sample_flags: sample identifier with values: 's', 'S', 'b', 'B',
          'k', 'K'.
       o  external_sample_name: the name of the sample supplied by the sample
          submitter.
       o  tissue: the name of the tissue from which the sample was collected
       o  wrap_group: the name of the group to which the results will
          be distributed. This information is used by the bbi-scirna-wrap
          script.
       o  lanes: the flowcell lanes that have the sample (library)
       o  hash_file: the path to a text file that has the hash barcode
          sequences for a sciPlex experiment.
       o  rt_file: path to custom RT primer barcode file.
       o  ligation_file: path to custom ligation barcode file.
       o  p7_file: path to custom P7 barcode file.
       o  p5_file: path to custom P5 barcode file.
       o  library: an arbitrary string, which can label the library or run.
          The intention is to store the string in the CDS.
       o  process_group: an integer that distinguishes sets of samples that
          are to be processed separately, in the sense that they result in
          different expression matrices, cdses, etc. For example,
          different libraries that are run in different Novaseq X lanes or
          have different sets of PCR primer pairs. Each sample in a process
          group is given the same process_group integer value. The samples in
          the first group are given the value 1 and samples in additional
          groups are given successive integer values. It is important to set
          the process_group values correctly when a run has more than one
          library where different libraries may have samples with the same
          names. It is critical to use distinct process_group values when
          the libraries in different lanes use different barcode sets; that
          is, different barcode files.
  o  sample names:
       o  begin with an alphanumeric character: a-z, A-Z, and 0-9
       o  allowed characters are alphabetic (a-z and A-Z), numeric (0-9), and
          period '.'
       o  this program converts other characters in the sample name to '.',
          and then checks for sample name degeneracy. If there is, the program
          exits immediately.
  o  external sample names:
       o  must not contain non-printable characters and no single or double
          quotes or backslashes. Cells may be empty.
  o  tissue:
       o  must not contain non-printable characters and no single or double
          quotes or backslashes. Cells may be empty.
  o  genomes:
       o  recognized genome names are listed in the variable 'genome_name_list'
          in this program's code. If a samplesheet genome name is not in the
          list, scirna_samplesheet.py gives a warning in case the name is
          mis-spelled.
       o  genome is the name of the organism that was sequenced. This
          identifies the files required to analyze the reads. The genome string
          is passed to the processing pipeline, and the pipeline uses it to
          find the required genome files. The file system paths to the
          available genomes are defined in the
          bbi-scirna-analyze/bin/star_genomes.txt file.
  o  wrap_group:
       o  valid names consist of lower and upper case alphabetic characters,
          numerals, and the symbols '.', '_', and '-'. Cells may be empty in
          which case those samples are excluded from the wrap output.
  o  wells:
       o  samplesheet wells are converted to indices where indices refer to
          physical wells, which are in the following order: rt and p7
          indices increase by column number along each row and p5 indices
          increase by row letter down each column. This is important when
          specifying wells by ranges. For example, the p7 'A' row is given
          by A01:A12. The p7 '1' column must be specified by each well; that
          is A01,B01,C01,D01,E01,F01,H01 (NOT A01:H01, which would give
          columns 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, and the first well in
          the 12th column.) Similary, the p5 '1' column is given by A01:H01
          and the p5 'A' row is necessarily given as A01,A02,A03,A04,A05,
          A06,A07,A08,A09,A10,A11,A12 (NOT A01:A12).
       o  the p7_well and p5_well or p7_rows and p5_columns is set to '0'
          (zero)  when the corresponding index is not detected by the
          sequencing machine. If one p7 well is set to '0', all p7 wells
          must be '0'. The same is true for p5 wells, p7 rows, and p5
          columns.
       o  wells that include a plate identifier must have '-' separating the
          plate and well,
          examples:
            o  P1-A10
            o  P02-H07
          The plate number consists of one or two digits.
       o  well ranges are given as a first and last well in a series of
          contiguous wells. The two well names are separated by a colon ':'.
          A row of rt or p7 wells is A01:A12 and a column of p5 wells is
          A01:H01.
          examples:
            o  A01:A12
            o  P01-B05:P01-B12
       o  wells are given as single or ranges of wells separated by commas.
          examples:
            o  A01:A12,B01:B12
            o  P1-A01,P1-A12
       o  well names do not depend on case
          examples:
            o  P1-A10
            o  p1-a10
  o  p7 rows:
       o  u-titer plate p7 rows used for PCR reactions
       o  rows are given as single rows separated by commas
          examples:
            o  E,F,G
  o  p5 columns:
       o  u-titer plate p5 columns used for PCR reactions
       o  columns are given as single columns separated
          by commas
          examples:
            o  3,4,5
            o  none
  o  p7 rows and p5 columns:
       o  by default, the i-th p7 row is paired with the i-th p5 column. No
          other combinations of specified rows and columns are considered.
          Consequently, the order of the rows and columns must be correct.
          Technical notes:
            o  in order to restrict the PCR index combinations to one p7 row
               and one p5 column per reaction plate, this program expands the
               number of spreadsheet rows per sample by the number of
               p7 rows/p5 columns pairs given for the sample in the input CSV
               file. For example, if for a sample, the p7_rows specification
               is 'A,B' and the p5_columns is '6,5', this program outputs two
               rows for the sample, the first with PCR row/column pair A/6
               and the second with PCR row/column pair B/5.
  o  p7 wells and p5 wells
       o  this program converts wells to indices and passes the indices to the
          pipeline. The pipeline uses combinations of all p7 and p5 wells. The
          order of the wells is unimportant. Use the string '0' to ignore the
          PCR index.
  o  lanes:
       o  given as one or more lanes or ranges of lanes using a colon, ':',
          to specify lane ranges, and commas, ',', to separate lanes or lane
          ranges. The lanes column may be omitted in which case all lanes will
          be used for all samples. If the lane column is included, all column
          cells must have values.
          examples:
            o  1:3
            o  1,3,5:8
  o  hash_file:
       o  path to the hash barcode file for samples with hash reads; leave the
          cell empty for samples without hash reads. Different samples in a run
          may have different hash_file values.
  o  sample_flags:
       o  sample_flags identify certain sample types, for example, Sentinel and
          Barnyard. The allowed values are 's', 'S', 'b', 'B', 'k', and 'K'.
          The sentinel samples are identified by an 's' the barnyard samples
          are identified by a 'b', and the keyhole experiment is identified by
          'k'. The samples with capitalized letters are used in the dashboard
          sentinel, barnyard, and keyhole reports. A column cell may have
          up to one flag character. Cells may be empty.
  o  rt_file:
       o  paths to custom rt barcode file. The path is not checked for
          validity. Empty cells denote use of the default barcode file.
       o  samples in a process_group must have the same rt_file.
       o  samples in a lane must have the same rt_file.
  o  ligation_file:
       o  paths to custom ligation barcode file. The path is not checked for
          validity. Empty cells denote use of the default barcode file.
       o  samples in a process_group must have the same ligation_file.
       o  samples in a lane must have the same ligation_file.
  o  p7_file:
       o  paths to custom p7 barcode file. The path is not checked for
          validity. Empty cells denote use of the default barcode file.
       o  samples in a process_group must have the same p7_file.
       o  samples in a lane must have the same p7_file.
  o  p5_file:
       o  paths to custom p5 barcode file. The path is not checked for
          validity. Empty cells denote use of the default barcode file.
       o  samples in a process_group must have the same p5_file.
       o  samples in a lane must have the same p5_file.
  o  library:
       o  library an arbitrary string that can used to identify the sequencing
          library
  o  process_group:
       o  a process_group is a single integer value
       o  each sample in a group is given the same integer value
       o  the samples in the first group are given the value 1
       o  samples in additional groups are given successive integer
          values.

  Notes:
    o  the well value sets are enclosed in quotes in this example because the
       sets may have commas.  Do not use quotes in the spreadsheet program
       cells; the spreadsheet program adds them when it writes the CSV file.
    o  the '-n' command line parameter identifies the number of flowcell
       lanes.
    o  when the lanes column is absent, it is assumed that all samples are
       present in all lanes.
    o  these pipelines store reads in BAM files (unaligned and aligned).
    o  the general work flow is
         o  bcl-convert makes fastq files demultiplexed by PCR barcode pairs.
            There is one pair of fastq files (fwd and rev reads) for each
            valid PCR barcode pair (and the fastq filenames have the p7 and p5
            barcode sequence indices, as well as the lane number, in them).
            These fastq files are used internally only, and are not returned
            to the user, i.e., 'published' by the Nextflow pipeline.
         o  bcl-convert samplesheet file identifies samples by lane, p7, and
            p5 sequences. It has no knowledge of process groups, rt barcodes,
            and ligation barcodes. Libraries that use different PCR oligo
            primer sets (or primer oligos in different wells), must be
            sequenced in different lanes.
         o  the reads in these fastq files are demultiplexed by RT and ligation
            barcodes. The RT barcodes identify reads by sample. The resulting
            reads are written to unaligned BAM files where all reads in a BAM
            file have the same lane and sample, p7, and p5 barcodes. The lane
            number and PCR barcode sequence indices are part of the BAM file
            name. The BAM file names have the form
              <sample_name>-<lane_number>_<p7_index>_<p5_index>-L<lane_index>.bam
            They are 'published' to the 'demux_out' directory.
         o  demultiplexed BAM files are merged by lanes, as required. This
            gathers reads by sample, PCR pair, and process_group. Their names
            have the form
              <sample_name>-<process_group>_<p7_index>_<p5_index>.merged.bam
            They are not published.
         o  if you have the same sample, with the same sample_name, in
            multiple libraries, each library run in a different lane(s), and
            you need the pipeline to process separately the reads from each
            library, assign a distinct process_group value to the sample
            rows for each library.
         o  the unaligned BAM files are processed by trimgalore to trim off
            adapter sequence. The resulting BAM file names have the form
              <sample_name>-<process_group>_<p7_index>_<p5_index>.trimmed.bam
            These BAM files are not published. 
         o  the trimmed read BAM files are aligned to the reference genome and
            assigned to cells by STARsolo. The STARsolo output files are for
            reads that have the same sample_name, process_group, and p7 and
            p5 barcodes. These aligned read BAM files are not 'published'.
         o  these (sample, process_group, pcr pair) aligned output files are
            merged such that all reads with the same sample_name and
            process_group are written to the same output BAM file. The merged
            BAM, count matrix, and statistics files are 'published' to the
            'analyze_out' directory.
            The BAM file names have the form
              <sample_name>-<process_group>.aligned.bam
         o  a cds file and umap.png file are made for each sample and
            process_group expression count matrix. These files are 'published'
            to the 'analyze_out' directory.
         o  when the hash_file value is set for a sample, the untrimmed BAM
            files are processed to find candidate hash reads and a cds is
            made with the hash read information.

Command line options:

usage: scirna_samplesheet.py [-h] [-i INPUT] [-o OUTPUT] [-s SEQUENCER_CLASS]
                             [-r RUN_DIR] [-l {2,3}] [-n NUMBER_LANES] [-e]
                             [-d] [-v]

A program to convert sci-RNA CSV samplesheet to pipeline samplesheet.

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input CSV samplesheet filename (required string).
  -o OUTPUT, --output OUTPUT
                        Output samplesheet filename (required string).
  -s SEQUENCER_CLASS, --sequencer_class SEQUENCER_CLASS
                        Sequencing machine class (not required string:
                        'illumina' or 'ultima'. Default: illumina
  -r RUN_DIR, --run_dir RUN_DIR
                        Illumina run directory path (optional string).
  -l {2,3}, --level {2,3}
                        Two or three level sci-rna-seq experiment (default:
                        3) (optional integer).
  -n NUMBER_LANES, --number_lanes NUMBER_LANES
                        Number of flowcell lanes (required integer).
  -e, --template        Write template samplesheet file
                        ('samplesheet.template.csv') with standard column
                        formats and exit (optional flag).
  -d, --documentation   Display documentation and exit (optional flag).
  -v, --version         Give program and JSON output file versions and exit
                        (optional flag).


  -- Example simple samplesheet file:

  Notes:
    o  all lanes have the same library (one process_group)
    o  required columns only

  sample_name,genome,rt_wells,p7_wells,p5_wells
  sample.1,mouse,"P01-A01:P01-A02","C01:C12","A11:H11"
  sample.2,mouse,"P01-A03:P01-A04","C01:C12","A11:H11"
  Barnyard,barnyard,"P01-A05:P01-A06","C01:C12","A11:H11"


  -- Example multi-library samplesheet file:

  Notes: 
    o  two libraries in two lane sets (two process_groups)
    o  the sample names in the two libraries are the same
    o  required columns only

  sample_name,genome,rt_wells,p7_wells,p5_wells,lanes,process_group
  my_sample,mouse,"P01-A01:P01-A02","C01:C12","A11:H11",1,1
  my_sample,mouse,"P01-A03:P01-A04","C01:C12","A11:H11",2,2
  Barnyard,barnyard,"P01-A05:P01-A06","C01:C12","A11:H11",1,1
  Barnyard,barnyard,"P01-A05:P01-A06","C01:C12","A11:H11",2,2

  -- Example complex multi-library samplesheet file with distinct
     lanes and PCR pairs to distinguish samples in process_groups

  Notes:
    o  three libraries in two lane sets (three process_groups)
    o  the sample names in the three libraries are the same
    o  required columns only

  sample_name,genome,rt_wells,p7_wells,p5_wells,lanes,process_group
  my_sample,mouse,"P01-A01:P01-A02","C01:C12","A11:H11",1,1
  my_sample,mouse,"P01-A03:P01-A04","C01:C12","A11:H11",2,2
  my_sample,mouse,"P01-A03:P01-A04","C01:C12","A12:H12",1:2,3
  Barnyard,barnyard,"P01-A05:P01-A06","C01:C12","A11:H11",1,1
  Barnyard,barnyard,"P01-A05:P01-A06","C01:C12","A11:H11",2,2
  Barnyard,barnyard,"P01-A05:P01-A06","C01:C12","A12:H12",1:2,3
"""

#
# Adding a new column:
#   o  add to column_header_name_list, columns_required_list or columns_optional_list, and, maybe column_allow_empty_cell
#   o  edit parse_header_column_name() to add condition in if elif ... statement to give value to column_name_dict
#   o  edit make_samplesheet_indexes() to add to row_out_list.append()  # add an elif() too
#   o  edit write_samplesheet_json_format() to add to sample_index_list.append() or sample_data dict (depending on where to store the values in the JSON file)
#   o  add default value if needed
#   o  add checks function
#


import sys
import re
import csv
import json
import argparse
import platform
import copy


#
# Samplesheet JSON file version.
#
program_version = '0.1.8'
json_file_version = '1.0.0'


#
# Default barcode file paths.
#
default_rt_file = ''
default_ligation_file = ''
default_p7_file = ''
default_p5_file = ''


#
# List of recognizable genome names.
# This program issues a warning if a samplesheet genome is
# not in this list.
#
genome_name_list = [
  'arabidopsis',
  'barnyard',
  'bat',
  'cat',
  'chicken',
  'corn',
  'cow',
  'cynomolgus',
  'dog',
  'drosophila',
  'duck',
  'elephant',
  'horse',
  'human',
  'macaque',
  'mouse',
  'opossum',
  'pig',
  'rabbit',
  'rat',
  'snake',
  'worm',
  'zebrafish',
  'hg19',
  'mm19',
  'hg19_mm9'
]


#
# Barcode file types.
#
barcode_file_types = [ 'rt_file', 'ligation_file', 'p7_file', 'p5_file' ]

#
# List of recognizable CSV column header names.
# These are used to check labels in the file.
#
p5_column_names = 'p5_wells p5_columns'
p7_column_names = 'p7_wells p7_rows'

column_header_name_list = [ 'rt_wells', 'sample_name', 'genome', 'sample_flags', 'external_sample_name', 'tissue', 'wrap_group', 'lanes', 'hash_file', 'rt_file', 'ligation_file', 'p7_file', 'p5_file', 'library', 'process_group' ]
column_header_name_list.extend( p5_column_names.split() )
column_header_name_list.extend( p7_column_names.split() )

columns_required_list = [ 'rt', 'p5', 'p7', 'sample_name', 'genome' ]
columns_optional_list = [ 'sample_flags', 'external_sample_name', 'tissue', 'wrap_group', 'lanes', 'hash_file', 'rt_file', 'ligation_file', 'p7_file', 'p5_file', 'library', 'process_group' ]

#
# Columns that may have empty cells. This is used to prevent
# errors when testing for empty cells in check_rows().
#
column_allow_empty_cell = [ 'sample_flags', 'external_sample_name', 'tissue', 'wrap_group', 'hash_file', 'rt_file', 'ligation_file', 'p7_file', 'p5_file', 'library' ]


#
# P7 row regex and P5 column regex.
#
p7_re_row_pattern_1 = r'([a-hA-H])([-:]([a-hA-H]))?$'
p5_re_col_pattern_1 = r'([1-9][0-2]?)([-:]([1-9][0-2]?))?$'
p7_re_row_pattern_2 = r'([a-hA-H0])([-:]([a-hA-H0]))?$'
p5_re_col_pattern_2 = r'([0-9][0-2]?)([-:]([0-9][0-2]?))?$'


#
# sample_flags regex.
#
sample_flags_re_pattern = r'([sSbBkK])$'


# OK
def display_documentation():
  print( __doc__ ) 
  return( 0 )
    

# OK
def clean_string(instring):
  # Remove leading and trailing whitespace.
  outstring = instring.strip()
  # Convert contiguous tabs to a space.
  outstring = re.sub( r'[\t]+', ' ', outstring )
  return( outstring )


# OK
def is_printable(string):
  # Check for non-printable characters.
  if(not string.isprintable()):
    return(False)
  return( True )


# OK
def has_space(string):
  if( re.search( r'[ ]', string ) ):
    return( True )
  return( False ) 


# OK
# Test == 1
#   Check that string has none of the following
#   characters:
#     backspace
#     form feed
#     newline
#     carriage return
#     single quote
#     double quote
#     backslash
# Test == 2
#   Check that string has only the following
#   characters:
#     alphabetic
#     numeric
#     ._-
def test_string(string, test):
  if( int( test ) == 1 ):
    if( not re.search( r'[\b\f\n\r\\\'\"]', string ) ):
      return( True )
  elif( int( test ) == 2 ):
    if( not re.search( r'[^a-zA-Z0-9._-]', string ) ):
      return(True)
  else:
    print( 'Error: test_string: unknown test value: %d' % ( int( test ) ) )
    return( False )
  return( False)


# OK
#
# Trim off leading and trailing whitespace and convert
# contiguous tabs to a space.
#
def clean_samplesheet_data( column_name_list, samplesheet_row_list ):
  errorFlag = False
  for irow, row_elements in enumerate(samplesheet_row_list):
    for icol in range( len( row_elements ) ):
      clean_string( row_elements[icol] )
      if( not is_printable( row_elements[icol] ) ):
        print('Error: non-printable character(s) in row %d column %d' % ( irow, icol + 1 ))
        errorFlag = True
  if( errorFlag ):
    sys.exit( -1 )


def check_args( args ):
  """
  There are no argument checks at this time.
  """
  error_string = ''
  if( len( error_string ) > 0 ):
    print( 'Command line argument errors:' )
    print( error_string, file=sys.stderr )
    sys.exit( -1 )
  return( 0 )


# OK
def parse_header_column_name( string_in, column_name_list, error_string ):
  """
  Split column header name into a 'type' and a 'format' and store as dictionary in column_name_list.
  """
  if( not string_in.lower() in column_header_name_list ):
    error_string += '  %s' % ( string_in )
    return( column_name_list, error_string )
  string_in = string_in.lower().strip()
  if( string_in == 'sample_name' ):
    column_name_dict = { 'type': 'sample_name', 'format': None }
  elif( string_in == 'genome' ):
    column_name_dict = { 'type': 'genome', 'format': None }
  elif( string_in == 'sample_flags' ):
    column_name_dict = { 'type': 'sample_flags', 'format': None }
  elif( string_in == 'external_sample_name' ):
    column_name_dict = { 'type': 'external_sample_name', 'format': None }
  elif( string_in == 'tissue' ):
    column_name_dict = { 'type': 'tissue', 'format': None }
  elif( string_in == 'wrap_group' ):
    column_name_dict = { 'type': 'wrap_group', 'format': None }
  elif( string_in == 'lanes' ):
    column_name_dict = { 'type': 'lanes', 'format': None }
  elif( string_in == 'hash_file' ):
    column_name_dict = { 'type': 'hash_file', 'format': None }
  elif( string_in == 'rt_file' ):
    column_name_dict = { 'type': 'rt_file', 'format': None }
  elif( string_in == 'ligation_file' ):
    column_name_dict = { 'type': 'ligation_file', 'format': None }
  elif( string_in == 'p7_file' ):
    column_name_dict = { 'type': 'p7_file', 'format': None }
  elif( string_in == 'p5_file' ):
    column_name_dict = { 'type': 'p5_file', 'format': None }
  elif( string_in == 'library' ):
    column_name_dict = { 'type': 'library', 'format': None }
  elif( string_in == 'process_group' ):
    column_name_dict = { 'type': 'process_group', 'format': None }
  else:
    mobj = re.match( r'([p][57]|rt)_(wells|rows|columns)', string_in )
    column_name_dict = { 'type': mobj.group( 1 ), 'format': mobj.group( 2 ) }
  column_name_list.append( column_name_dict )
  return( column_name_list, error_string )


# OK
def check_header_column_names( column_name_list ):
  """
  Check column header names.
  Notes:
    o  check that required columns occur
    o  check that each allowed column type occurs once
  """
  columns_allowed = {}
  for column_name in columns_required_list:
    columns_allowed[column_name] = 0
  for column_name in columns_optional_list:
    columns_allowed[column_name] = 0
  for column_name_dict in column_name_list:
    columns_allowed[column_name_dict['type']] += 1
  errorFlag = 0
  # Check for rt, p5, and p7 specification columns.
  for column_name in columns_allowed:
    if( columns_allowed[column_name] == 0 and column_name not in columns_optional_list ):
      print( 'Error: column for \'%s\' missing.' % ( column_name ), file=sys.stderr )
      if( column_name == 'rt' ):
        print( '  acceptable rt header values: %s' % ( rt_column_names ), file=sys.stderr )
      elif( column_name == 'p5' ):
        print( '  acceptable p5 header values: %s' % ( p5_column_names ), file=sys.stderr )
      elif( column_name == 'p7' ):
        print( '  acceptable p7 header values: %s' % ( p7_column_names ), file=sys.stderr )
      errorFlag = 1
    elif( columns_allowed[column_name] > 1 ):
      print( 'Error: column for \'%s\' occurs %d times.' % ( column_name, columns_allowed[column_name] ), file=sys.stderr )
      errorFlag = 1
  if( errorFlag ):
    sys.exit( -1 )
  p5_col = False
  p7_row = False
  for column_name_dict in column_name_list:
    if( column_name_dict['type'] == 'p5' and column_name_dict['format'] == 'columns' ):
      p5_col = True
    if( column_name_dict['type'] == 'p7' and column_name_dict['format'] == 'rows' ):
      p7_row = True
  if( p5_col != p7_row ):
    print( 'Error: p5 is %sin \'columns\' format but p7 is %sin \'rows\' format' % ( '' if p5_col else 'not ', '' if p7_row else 'not ' ), file=sys.stderr )
    sys.exit( -1 )
  return( 0 )


# OK
def parse_header( row_header ):
  """
  Convert column header (row) into a list of column name dictionaries.
  The dictionary has the elements
    key     value description
    type    entry type name: rt, p5, p7, sample_name, genome, sample_flags, external_sample_name, tissue, wrap_group, lanes, hash_file
    format  barcode format values: wells, None (see column_header_name_list for allowed combinations of type and format)
  """
  column_name_list = []
  error_string = ''
  for str in row_header:
    column_name_list, error_string = parse_header_column_name( str, column_name_list, error_string )
  if(len(error_string) > 0):
    print('Error: invalid header label(s): %s' % (error_string))
    sys.exit(-1)
  check_header_column_names( column_name_list )
  return( column_name_list )


# OK
def well_to_index( plate, row, column, across_row_first=True, element_coordinates=[None,None] ):
  """
  Convert a well specification to a plate index.  max was 384

  Args:
    plate              integer plate number >= 1
    row                character row (A-H)
    column             integer column number (1-12)
    across_row_first   bool index increases by one as a row is traversed; that is,
                       moving from column to column along row
  Returns:
    index: an integer well index

  """
  row = row.lower()
  if( plate < 1 or
      not row in [ 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h' ] or
      column < 1 or column > 12 ):
    print( 'Error: spreadsheet cell: %s%s:  bad well values: plate: %d  row: %s  col: %d' % ( element_coordinates[0], element_coordinates[1], plate, row, column ), file=sys.stderr )
    sys.exit( -1 )
  irow = ord( row ) - ord( 'a' )
  icol = column - 1
  if( across_row_first ):
    well_index = irow * 12 + icol + 1
  else:
    well_index = icol * 8 + irow + 1
  return( well_index + ( plate - 1 ) * 96 )


# OK
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


# OK
def check_index_list( index_type, index_list, element_coordinates = [ None, None ] ):
  """
  Check and clean index list
    o  check for duplicate indices
    o  remove duplicate indices
  """
  index_dict = {}
  for i in index_list:
    index_dict.setdefault( i, 0 )
    index_dict[i] += 1
  duplicate_list = []
  for i in index_dict.keys():
    if( index_dict[i] > 1 ):
      duplicate_list.append( str( i ) )
  if( len( duplicate_list ) > 0 ): 
    print( 'Warning: spreadsheet cell: %s %s: duplicate %s index(es): %s' % ( index_type, element_coordinates[0], element_coordinates[1], ' '.join( sorted(duplicate_list) ) ), file=sys.stderr )
    print( '         These indices are duplicated within the reported spreadsheet cell, which may be intentional.')
  return( list( set( index_list ) ) ) 


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


# OK
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


# OK probably
def parse_wells( string_in, across_row_first=True, element_coordinates = [ None, None ] ):
  """
  Convert a well specification to an index string.
  Acceptable well specifications include
    o  a single well without a plate specified (implicit plate=1), examples
         o  A5
         o  A05
    o  a single well with a plate specified
         o  P1-A5
    o  range of wells without plates specified,
       Note: the range of indices depends on whether the reaction type is n5/p5 or n7/p7.
         o  A9:B9
    o  range of wells with plates specified,
         o  P1-A5:P1-A12
    o  single and/or ranges of wells separated by commas
         o  A10,P1-B5:P1-B10,C7
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for well_range in string_in.split( ',' ):
#    mobj =  re.match( r'([pP][0]?([1-4*])[-])?([a-hA-H])([0]?[1-9][0-2]?)([:]([pP][0]?([1-4*])[-])?([a-hA-H])([0]?[1-9][0-2]?))?$', well_range )
    if( mobj :=  re.match( r'([pP][0]?([1-8])[-])?([a-hA-H])([0]?[1-9][0-2]?)([:]([pP][0]?([1-8])[-])?([a-hA-H])([0]?[1-9][0-2]?))?$', well_range.strip() ) ):
      # first well
      row1 = mobj.group( 3 )
      col1 = int( mobj.group( 4 ) )
      if( col1 < 1 or col1 > 12 ):
        print( 'Error: spreadsheet cell: %s%s: bad well: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
      # is plate specified?
      if( mobj.group( 1 ) ):
        plate1_list = [ int( mobj.group( 2 ) ) ]
      else:
        plate1_list = [ 1 ]
      #
      if( mobj.group( 5 ) ):
        if( ( mobj.group( 2 ) == None ) != ( mobj.group( 7 ) == None ) ):
          print( 'Error: spreadsheet cell: %s%s: either both or neither well in a range must have plates specified: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
          sys.exit( -1 )
        # second well, if this is a range
        row2 = mobj.group( 8 )
        col2 = int( mobj.group( 9 ) )
        if( col2 < 1 or col2 > 12 ):
          print( 'Error: spreadsheet cell: %s%s: bad well: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
          sys.exit( -1 )
        # is plate specified?
        if( mobj.group( 6 ) ):
          plate2_list = [ int( mobj.group( 7 ) ) ]
        else:
          plate2_list = [ 1 ]
        #
      else:
        plate2_list = plate1_list
        row2 = row1
        col2 = col1
      #
      for plate1, plate2 in zip( plate1_list, plate2_list ):
        index1 = well_to_index( plate1, row1, col1, across_row_first, element_coordinates )
        index2 = well_to_index( plate2, row2, col2, across_row_first, element_coordinates )
        if( index2 < index1 ):
          print( 'Error: spreadsheet cell: %s%s: bad well range: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
          sys.exit( -1 )
        for i in range( index1, index2 + 1 ):
          index_list.append( i )

    elif( mobj := re.match( r'0$', well_range.strip() ) ):
      index_list.append( 0 )
    else:
      if( not mobj ):
        print( 'Error: spreadsheet cell: %s%s: bad well or well range \'%s\'' % ( element_coordinates[0], element_coordinates[1], well_range ), file=sys.stderr )
        sys.exit( -1 )
  return( check_index_list( 'well', index_list, element_coordinates ) )


def parse_rows( string_in, element_coordinates = [ None, None ] ):
  """
  Convert a row specification to an index string.
  Acceptable row specifications include
    o  single row
         o  B
    o  row range
         o  E-G
    o  single and/or ranges of rows separated by commas
         o  E-F,H
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for row_range in string_in.split( ',' ):
    if( mobj := re.match( p7_re_row_pattern_1, row_range.strip() ) ):
      row1 = mobj.group( 1 )
      row1_index = well_to_index( 1, row1, 1, True, element_coordinates )
      row2_index = row1_index
      if( mobj.group( 2 ) ):
        row2 = mobj.group( 3 )
        row2_index = well_to_index( 1, row2, 1, True, element_coordinates )
        if( row2_index < row1_index ):
          print( 'Error: spreadsheet cell: %s%s: bad row range: \'%\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
          sys.exit( -1 )
      index1 = row1_index
      index2 = row2_index + 11
      for i in range( index1, index2 + 1 ):
        index_list.append( i )
    elif( mobj := re.match( r'0$', row_range.strip() ) ):
      # '0' means that the barcode was not sequenced
      index_list.append( 0 )
    else:
      if( not mobj ):
        print( 'Error: spreadsheet cell: %s%s: bad row or row range \'%s\'' % ( element_coordinates[0], element_coordinates[1], row_range ), file=sys.stderr )
        sys.exit( -1 )
  return( check_index_list( 'row', index_list, element_coordinates ) )


def parse_columns( string_in, element_coordinates = [ None, None ] ):
  """
  Convert a column specification to an index string.
  Acceptable column specifications include
    o  single column
         o  5
    o  column range
         o  6-8
    o  single and/or ranges of column separated by commas
         o  9-11,3
    o  column values must be integers in the range 0 to 12
  """
  index_list = []
  string_in = re.sub( r'\s', '', string_in )
  for col_range in string_in.split( ',' ):
    if(col_range == 'none'):
      index_list.append(0)
      continue
    if( mobj := re.match( p5_re_col_pattern_1, col_range.strip() ) ):
      col1 = int( mobj.group( 1 ) )
      col1_index = well_to_index( 1, 'A', col1, False, element_coordinates )
      col2_index = col1_index
      if( col1 < 1 or col1 > 12 ):
        print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
        sys.exit( -1 )
      if( mobj.group( 2 ) ):
        col2 = int( mobj.group( 3 ) )
        if( col2 < 1 or col2 > 12 ):
          print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
          sys.exit( -1 )
        if( col2 < col1 ):
          print( 'Error: spreadsheet cell: %s%s: bad column range: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
          sys.exit( -1 )
        col2_index = well_to_index( 1, 'A', col2, False, element_coordinates )
      index1 = col1_index
      index2 = col2_index + 7
      for i in range( index1, index2 + 1 ):
        index_list.append( i )
    elif( mobj := re.match( r'0$', col_range.strip() ) ):
      # '0' means that the barcode was not sequenced
      index_list.append( 0 )
    else:
      if( not mobj ):
        print( '2 Error 1: spreadsheet cell: %s%s: bad column or column range \'%s\'' % ( element_coordinates[0], element_coordinates[1], col_range ), file=sys.stderr )
        sys.exit( -1 )
  return( check_index_list( 'column', index_list, element_coordinates ) )


def parse_sample_flags(string_in, element_coordinates ):
  sample_flag_list = []
  for sample_flag in string_in.split( ',' ):
    if( len( sample_flag ) == 0 ):
      return( [ '' ] )
    mobj = re.match( sample_flags_re_pattern, sample_flag.strip() )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad sample flag specification: \'%s\'' % ( element_coordinates[0], element_coordinates[1], sample_flag ) )
      sys.exit( -1 )
    sample_flag_list.append( sample_flag )
  return( sample_flag_list )


# lane specification
#   o  lanes are integers values > 0
#   o  there are one or more lane and lane_ranges
#   o  lanes and lane ranges are separated by a comma
#
def parse_lanes(string_in, element_coordinates ):
  lane_list = []
  for lane_range in string_in.split( ',' ):
    mobj = re.match(r'([0-9]?)([-:]([0-9]))?$', lane_range.strip())
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad lane specification: \'%s\'' % ( element_coordinates[0], element_coordinates[1], lane_range ) )
      sys.exit( -1 )
    if( mobj.group( 1 ) != '' ):
      lane1_index = int( mobj.group( 1 ) )
    else:
      return( '' )
    lane2_index = lane1_index
    if( lane1_index < 0 ):
      print( 'Error: spreadsheet cell: %s%s: bad lane value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], lane1_index ), file=sys.stderr)
      sys.exit( -1 )
    if( mobj.group( 2 ) ):
      lane2_index = int( mobj.group( 3 ) )
      if( lane2_index <= lane1_index ):
        print( 'Error: spreadsheet cell: %s%s: bad lane value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], lane2_index ), file=sys.stderr)
        sys.exit( -1 )
    for i in range( lane1_index, lane2_index + 1 ):
      lane_list.append( i )
  return( check_index_list( 'lane', lane_list, element_coordinates ) )


# OK
def check_rows( column_name_list, csv_rows ): 
  """
  Trim off empty cells at end of rows and columns
  and check for internal empty cells. Allow empty
  internal row.
  Notes:
    o  we expect
         o  nrows
  """
  # check for internal empty row 
  csv_rows_out = []
  num_col = len( column_name_list )
  for irow, row_elements in enumerate( csv_rows ):
    if(len(row_elements) < num_col):
      print('Error: missing %d cell(s) in row %d' % ( num_col - len(row_elements), irow))
      sys.exit( -1 )
    # Number of cells that are supposed to have
    # content but don't.
    num_empty = 0 
    row_elements_out = []
    for icol, cell in enumerate( row_elements ):
      if( icol == num_col ):
        break
      # Add valid cell to output row.
      if( len( cell ) > 0 or ( column_name_list[icol]['type'] in column_allow_empty_cell ) ):
        row_elements_out.append( cell )
      else:
        print( 'Error: unexpected empty cell in column \'%s\'' % (column_name_list[icol]['type']))
        num_empty += 1 
    # Save the row if there are no invalid empty cells.
    if( num_empty == 0 ):
      csv_rows_out.append( row_elements_out )
    elif( num_empty > 0 and num_empty < num_col ):
      # Note: we allow for an empty row with
      #       num_empty == num_col.
      print( 'Error: row %d has empty cells' % ( irow + 1 ) )
      sys.exit( -1 )
  return( csv_rows_out )


# OK
def read_samplesheet( file ):
  """
  Read CSV samplesheet input file.
  Notes:
    o  the first row in the file must have column header names.
    o  the column header names must be in the list 'column_header_name_list'.
    o  the column order is arbitrary.
  """
  samplesheet_row_list = [] 
  csv_rows = csv.reader( file, delimiter=',', quotechar='"')
  csv_rows = list( csv_rows )
  row_header = csv_rows[0]
  column_name_list = parse_header( row_header )
  csv_rows = check_rows( column_name_list, csv_rows )
  for row_elements in csv_rows[1:]:
    samplesheet_row_list.append( row_elements )
  return( column_name_list, samplesheet_row_list )


# OK probably
def check_sample_names( column_name_list, samplesheet_row_list ):
  """
  Check for name degeneracy.
  Check sample names for unacceptable characters and, if present, convert them to '.'.
  Sample names must begin with [a-zA-Z].
  Unacceptable characters are characters that are not a-z, A-Z, 0-9, and '.'
  Check for name degeneracy after substitutions.
  Check that the barnyard sample is labeled 'Barnyard'.
  """
  sample_name_in_dict = {}
  sample_name_out_dict = {}
  num_sample_name = 0
  for row_elements in samplesheet_row_list:
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i] 
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'sample_name' ):
        continue
      sample_name_in_dict.setdefault( element_string, 0 )
      sample_name_in_dict[element_string] += 1
      num_sample_name += 1
      mobj = re.match( r'[a-zA-Z0-9]', element_string.strip() )
      if( not mobj ):
        print( 'Error: sample names must begin with an alphanumeric character', file=sys.stderr )
        sys.exit( -1 )
      # Convert invalid characters to periods.
      row_elements[i] = re.sub( r'[^a-zA-Z0-9.]', '.', element_string )
#      sample_name_out_dict.setdefault( element_string, True )
      sample_name_out_dict.setdefault( row_elements[i], True )
  errorFlag = False
  for sample_name in sample_name_in_dict.keys():
    if( sample_name_in_dict[sample_name] > 1 ):
      print( 'Warning: sample name \'%s\' not unique. It is used %d times.' % ( sample_name, sample_name_in_dict[sample_name] ), file=sys.stderr )
  if( len( sample_name_out_dict ) != len( sample_name_in_dict ) ):
    print( 'Error: unacceptable names are not distinct after editing', file=sys.stderr )
    errorFlag = True
  if( errorFlag ):
    sys.exit( -1 )
  # Check barnyard sample name. (This is unnecessary, I believe.)
#   for row_elements in samplesheet_row_list:
#     for i in range( len( row_elements ) ):
#       column_name_dict = column_name_list[i]
#       element_string = row_elements[i]
#       if( column_name_dict['type'] != 'sample_name' ):
#         continue
#       mobj = re.search( r'barn', element_string.lower() )
#       if( mobj and element_string != 'Barnyard' ):
#         print( '**' )
#         print( '** Warning: barnyard sample name (%s) not \'Barnyard\'.' % ( element_string ), file=sys.stderr )
#         print( '**          Consider re-naming it to \'Barnyard\' for compatibility' )
#         print( '**          with the experiment dashboard.' )
#         print( '**' )
#       break
  return( samplesheet_row_list )


# OK
def check_genome_names( column_name_list, samplesheet_row_list ):
  """
  Check genome names and warn if not in our list.
  """
  bad_row_dict = {}
  errorFlag = False
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'genome' ):
        continue
      if( has_space(element_string) ):
        print( 'Error: genome name \'%s\' in row %d has one or more spaces' % ( element_string, irow ) )
        errorFlag = 1
      if( not row_elements[i] in genome_name_list ):
        bad_row_dict.setdefault( row_elements[i], True )
  if( len( bad_row_dict.keys() ) > 0 ):
    print( 'The following genomes are not in my list of known genomes (they may be mis-spelled or not in my list).', file=sys.stderr )
    for missing_genome_name in bad_row_dict.keys():
      print( '  \'%s\'' % ( missing_genome_name ), file=sys.stderr )
  if( errorFlag ):
    sys.exit( -1 )
  return( 0 )


# OK
def check_external_sample_name( column_name_list, samplesheet_row_list ):
  """
  Check that each external_sample_name is valid.
  """
  errorFlag = False
  sample_name_in_dict = {}
  sample_name_out_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'external_sample_name' ):
        continue
      # Allow empty external sample name cell.
      if( len( element_string ) == 0 ):
        continue
      sample_name_in_dict.setdefault( element_string, 0 )
      sample_name_in_dict[element_string] += 1
      if(not test_string( element_string, 1 ) ):
        print( 'Error: bad external sample name in row %d' % ( irow + 2 ), file=sys.stderr )
        errorFlag = True
      sample_name_out_dict.setdefault( element_string, True )
  for sample_name in sample_name_in_dict.keys():
    if( sample_name_in_dict[sample_name] > 1 ):
      print( 'Warning: external sample name \'%s\' not unique. It is used %d times.' % ( sample_name, sample_name_in_dict[sample_name] ), file=sys.stderr )
  if( len( sample_name_out_dict ) != len( sample_name_in_dict ) ):
    print( 'Error: unacceptable names are not distinct after editing', file=sys.stderr )
    errorFlag = True
  if( errorFlag ):
    sys.exit( -1 )
  return( samplesheet_row_list )


# OK
def check_tissue( column_name_list, samplesheet_row_list ):
  """
  Check that sample tissue strings are valid.
  """
  bad_row_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for icol in range( len( row_elements ) ):
      column_name_dict = column_name_list[icol]
      element_string = row_elements[icol]
      if( column_name_dict['type'] != 'tissue' ):
        continue
      # Allow empty external sample name cell.
      if( len( element_string ) == 0 ):
        continue
      if(not test_string(element_string, 1)):
        bad_row_dict.setdefault( irow + 2, element_string )
  if( len( bad_row_dict.keys() ) > 0 ):
    print('Error: invalid tissue names:')
    print('  Row\tCell')
    for irow in bad_row_dict.keys():
      print( '  %s\t\'%s\'' % ( irow, bad_row_dict[irow] ) )
    sys.exit( -1 )
  return( samplesheet_row_list )


# OK
def check_wrap_group( column_name_list, samplesheet_row_list ):
  """
  Check that wrap_group values are valid.
  """
  bad_row_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'wrap_group' ):
        continue
      if( len( element_string ) == 0 ):
        continue
      if( not test_string( element_string, 2 ) ):
        bad_row_dict.setdefault( irow + 2, element_string )
  if( len( bad_row_dict.keys() ) > 0 ):
    print('Sample wrap_group values have invalid characters.' )
    print('  Row\tCell')
    for irow in bad_row_dict.keys():
      print( '  %s\t\'%s\'' % ( irow, bad_row_dict[irow] ) )
    sys.exit( -1 )
  return( 0 )


def check_lanes( column_name_list, samplesheet_row_list ):
  """ 
  Check that lanes values are valid.
  """
  bad_row_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'lanes' ):
        continue
      if( len( element_string ) == 0 ):
        bad_row_dict.setdefault( irow + 2, element_string)
        continue
      for lane_range in element_string.split( ',' ):
        mobj = re.match( r'^([0-9]+)([:]([0-9]+))?$', lane_range.strip())
        if( mobj == None or ( mobj.group( 3 ) != None and int( mobj.group( 3 ) ) <= int( mobj.group( 1 ) ) ) ):
          bad_row_dict.setdefault( irow + 2, element_string )
  if( len( bad_row_dict.keys() ) > 0 ):
    print('Sample lanes values have invalid characters or are empty.' )
    print('  Row\tCell')
    for irow in bad_row_dict.keys():
      print( '  %s\t\'%s\'' % ( irow, bad_row_dict[irow] ) )
    sys.exit( -1 )
  return( 0 )


def check_hash_file( column_name_list, samplesheet_row_list ):
  """ 
  Check that hash_file values are valid.
  """
  bad_row_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'hash_file' ):
        continue
      if( len( element_string ) == 0 ):
        continue
      if( not check_path_validity( element_string ) ):
        bad_row_dict.setdefault( irow + 2, element_string )
  if( len( bad_row_dict.keys() ) > 0 ):
    print('Sample hash reads values have invalid characters.' )
    print('  Row\tCell')
    for irow in bad_row_dict.keys():
      print( '  %s\t\'%s\'' % ( irow, bad_row_dict[irow] ) )
    sys.exit( -1 )
  return( 0 )


def check_path_validity( path ):
  """
  Perform simple test for valid path string.
  """
  # Accept any string and allow file open check whether the file exists.
  return( True )


def check_barcode_files( column_name_list, samplesheet_row_list ):
  """
  Check that barcode file path strings are (likely to be) allowed.
  This may not be a comprehensive test, or a test at all.
  """
  bad_row_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( not column_name_dict['type'].lower() in barcode_file_types ):
        continue
      if( len( element_string ) == 0 ):
        continue
      if( check_path_validity( element_string ) == False ):
        bad_row_dict.setdefault( irow + 2, element_string )
  if( len( bad_row_dict.keys() ) > 0 ):
    print('Sample file path strings have invalid characters.' )
    print('  Row\tCell')
    for irow in bad_row_dict.keys():
      print( '  %s\t\'%s\'' % ( irow, bad_row_dict[irow] ) )
    sys.exit( -1 )
  return( 0 )


#
# Check lane- and sample-oriented value consistency.
# Notes:
#   o  use a multilayer dictionary that looks like
#
#  counter[<lane_set value>]['sample_name'][<sample_name>] = <count>
#  counter[<lane_set value>]['rt_file'][<rt_file>] = <count>
#  counter[<lane_set value>]['ligation_file'][<ligation_file>] = <count>
#  counter[<lane_set value>]['p7_file'][<p7_file>] = <count>
#  counter[<lane_set value>]['p5_file'][<p5_file>] = <count>
#  counter[<lane_set value>]['library'][<library>] = <count>
#
#  counter[<lane_set value>]['sample_name'][<sample_name>]['genome'][<genome>] = <count>
#  counter[<lane_set value>]['sample_name'][<sample_name>]['sample_flags'][<sample_flags>] = <count>
#  counter[<lane_set value>]['sample_name'][<sample_name>]['tissue'][<tissue>] = <count>
#  counter[<lane_set value>]['sample_name'][<sample_name>]['wrap_group'][<wrap_group>] = <count>
#  counter[<lane_set value>]['sample_name'][<sample_name>]['hash_file'][<hash_file>] = <count>
#
def counter_set_value(counter_dict, row_value_dict, column_type, counter_flag):
  if(row_value_dict.get(column_type, None) == None):
    return(counter_dict)
  value = counter_dict.setdefault(column_type, None)
  if(counter_flag == False):
    if(value == None):
      counter_dict[column_type] = dict()
    value = counter_dict[column_type].setdefault(row_value_dict[column_type], None)
    if(value == None):
      counter_dict[column_type][row_value_dict[column_type]] = dict()
    return(counter_dict)
  elif( counter_flag == True):
    if(value == None):
      counter_dict[column_type] = dict()
    value = counter_dict[column_type].setdefault(row_value_dict[column_type], 0)
    counter_dict[column_type][row_value_dict[column_type]] += 1
    return(counter_dict)

 
def check_lane_sample_consistency( column_name_list, samplesheet_row_list ):
  """
  Check consistency of attributes by lane and sample.
  Expect that the values for the following attributes are
  consistent for a lane set:
    o  rt_file
    o  ligation_file
    o  p7_file
    o  p5_file
    o  library

  Expect that the values for the following attributes are
  consistent for a combination of lane set and sample name:
    o  genome
    o  flags
    o  tissue
    o  wrap_group
    o  hash_file
  """
  # Initialize counters.
  counter = dict()
  # Count lane sets.
  for irow, row_elements in enumerate(samplesheet_row_list):
    row_value_dict = dict()
    for icol in range( len( row_elements ) ):
      column_name_dict = column_name_list[icol]
      column_type = column_name_list[icol]['type']
      element_string = row_elements[icol]
      row_value = row_value_dict.setdefault(column_type, element_string)

#    lanes_value = row_value_dict['lanes']
    lanes_value = row_value_dict.get('lanes', 'x')
    sample_name_value = row_value_dict['sample_name']

    value = counter.setdefault(lanes_value, None)
    if(value == None):
      counter[lanes_value] = dict()

    counter[lanes_value] = counter_set_value(counter[lanes_value], row_value_dict, 'sample_name', False)
    counter[lanes_value] = counter_set_value(counter[lanes_value], row_value_dict, 'rt_file', True)
    counter[lanes_value] = counter_set_value(counter[lanes_value], row_value_dict, 'ligation_file', True)
    counter[lanes_value] = counter_set_value(counter[lanes_value], row_value_dict, 'p7_file', True)
    counter[lanes_value] = counter_set_value(counter[lanes_value], row_value_dict, 'p5_file', True)
    counter[lanes_value] = counter_set_value(counter[lanes_value], row_value_dict, 'library', True)

    counter[lanes_value]['sample_name'][sample_name_value] = counter_set_value(counter[lanes_value]['sample_name'][sample_name_value], row_value_dict, 'genome', True)
    counter[lanes_value]['sample_name'][sample_name_value] = counter_set_value(counter[lanes_value]['sample_name'][sample_name_value], row_value_dict, 'sample_flags', True)
    counter[lanes_value]['sample_name'][sample_name_value] = counter_set_value(counter[lanes_value]['sample_name'][sample_name_value], row_value_dict, 'tissue', True)
    counter[lanes_value]['sample_name'][sample_name_value] = counter_set_value(counter[lanes_value]['sample_name'][sample_name_value], row_value_dict, 'wrap_group', True)
    counter[lanes_value]['sample_name'][sample_name_value] = counter_set_value(counter[lanes_value]['sample_name'][sample_name_value], row_value_dict, 'hash_file', True)

  # print(json.dumps(counter, indent=2))

  # Check counts.
  # Lane tests.
  print('  Check lane-oriented value consistency...')
  error_flag_1 = False
  lanes_test_list = ['p7_file', 'p5_file']
  for lanes in counter.keys():
    for lanes_test in lanes_test_list:
      if(not lanes_test in row_value_dict.keys()):
        continue
      if(len(counter[lanes][lanes_test].keys()) > 1):
        if(error_flag_1 == False):
          print('    more than one value for')
        print('      \'%s\' in lanes %s' % (lanes_test, lanes))
        error_flag_1 = True
  if(error_flag_1 == False):
    print('    ** no problems noticed **')

  # Lane + sample tests.
  print('  Check lane- and sample-oriented value consistency...')
  error_flag_2 = False
  sample_test_list = ['genome', 'sample_flags', 'tissue', 'wrap_group', 'hash_file']
  for lanes in counter.keys():
    for sample_test in sample_test_list:
      if(not sample_test in row_value_dict.keys()):
        continue
      for sample_name in counter[lanes]['sample_name'].keys():
        if(len(counter[lanes]['sample_name'][sample_name][sample_test]) > 1):
          if(error_flag_2 == False):
            print('    more than one value for')
          print('      \'%s\' for sample %s in lanes %s' % (sample_test, sample_name, lanes))
          error_flag_2 = True
  if(error_flag_2 == False):
    print('    ** no problems noticed **')

  if(error_flag_1 == True or error_flag_2 == True):
    sys.exit(0)


def check_process_groups(column_name_list, samplesheet_row_list ):
  # check whether there is a process_group column
  # if no, return
  process_group_flag = False
  for column_name in column_name_list:
    if(column_name['type'] == 'process_group'):
      process_group_flag = True
      break
  if(process_group_flag == False):
    return(0)

  # check that each cell has an integer value
  print('check process_groups: check for integer values') # bge
  bad_row_dict = {}
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'process_group' ):
        continue
      if( not re.match(r'^[0-9]+$', element_string)):
        bad_row_dict.setdefault( irow + 2, element_string )
        continue
  if( len( bad_row_dict.keys() ) > 0 ):
    print('Sample process_group values have invalid characters.' )
    print('  Row\tCell')
    for irow in bad_row_dict.keys():
      print( '  %s\t\'%s\'' % ( irow, bad_row_dict[irow] ) )
    sys.exit( -1 )

  # check that the integer values start with 1 and are
  # sequential
  process_group_set = set()
  for irow, row_elements in enumerate(samplesheet_row_list):
    for i in range( len( row_elements ) ):
      column_name_dict = column_name_list[i]
      element_string = row_elements[i]
      if( column_name_dict['type'] != 'process_group' ):
        continue
      process_group_set.add(int(element_string))
  process_group_list_sorted = sorted(list(process_group_set))
  if(min(process_group_list_sorted) != 1):
    print('Minimum process_group value is not 1')
    sys.exit(-1)
  if(process_group_list_sorted != list(range(min(process_group_list_sorted), max(process_group_list_sorted)+1))):
    print('Process group values are not sequential')
    sys.exit(-1)
  return(0)


def expand_rows( string_in, element_coordinates = [ None, None ] ):
  """
  Expand a P7 row specification to a list of rows.
  Acceptable row specifications include
    o  single row
         o  B
    o  row range
         o  E-G
    o  single and/or ranges of rows separated by commas
         o  E-F,H
  """
  string_in = re.sub( r'\s', '', string_in )
  row_list = []
  for row_range in string_in.split( ',' ):
    mobj = re.match( p7_re_row_pattern_2, row_range.strip() )
    if( not mobj ):
      print( 'Error: spreadsheet cell: %s%s: bad row or row range \'%s\'' % ( element_coordinates[0], element_coordinates[1], row_range ), file=sys.stderr )
      sys.exit( -1 )
    icode_row1 = ord(mobj.group( 1 ))
    icode_row2 = icode_row1
    if( mobj.group( 2 ) ):
      icode_row2 = ord(mobj.group( 3 ))
      if( icode_row2 < icode_row1 ):
        print( 'Error: spreadsheet cell: %s%s: bad row range: \'%\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
    for icode in range(icode_row1, icode_row2+1):
      row_list.append(chr(icode))
  return(row_list)


def expand_columns( string_in, element_coordinates = [ None, None ] ):
  """
  Expand a P5 column specification to a list of columns.
  Acceptable column specifications include
    o  single column
         o  5
    o  column range
         o  6-8
    o  single and/or ranges of column separated by commas
         o  9-11,3
  """
  string_in = re.sub( r'\s', '', string_in )
  column_list = []
  for col_range in string_in.split( ',' ):
    if(col_range.strip() == 'none'):
      column_list.append('none')
      continue
    mobj = re.match( p5_re_col_pattern_2, col_range.strip() )
    if( not mobj ):
      print( '1 Error: spreadsheet cell: %s%s: bad column or column range \'%s\'' % ( element_coordinates[0], element_coordinates[1], col_range ), file=sys.stderr )
      sys.exit( -1 )

    icol_col1 = int( mobj.group( 1 ) )
    icol_col2 = icol_col1
    if( icol_col1 < 0 or icol_col1 > 12 ):
      print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
      sys.exit( -1 )
    if( mobj.group( 2 ) ):
      icol_col2 = int( mobj.group( 3 ) )
      if( icol_col2 < 0 or icol_col2 > 12 ):
        print( 'Error: spreadsheet cell: %s%s: bad column value: \'%d\'' % ( element_coordinates[0], element_coordinates[1], col1 ), file=sys.stderr )
        sys.exit( -1 )
      if( icol_col2 < icol_col1 ):
        print( 'Error: spreadsheet cell: %s%s: bad column range: \'%s\'' % ( element_coordinates[0], element_coordinates[1], string_in ), file=sys.stderr )
        sys.exit( -1 )
    for i in range( icol_col1, icol_col2 + 1 ):
      column_list.append(str(i))
  return(column_list)


def test_pcr_format(column_name_list):
  """
  Find the P7 and P5 columns in the header column_name_list.
  Return the column indices, or None if the P7 and P5 are
  not in the row and column format.
  """
  pcr7_column = None
  pcr5_column = None
  for icol in range(len(column_name_list)):
    column_dict = column_name_list[icol]
    if(column_dict['type'] == 'p7' and
       column_dict['format'] == 'rows'):
      pcr7_column = icol
    elif(column_dict['type'] == 'p5' and
         column_dict['format'] == 'columns'):
      pcr5_column = icol
  return(pcr7_column, pcr5_column)


def expand_sample_rows(column_name_list, samplesheet_row_list):
  """
  Expand sample rows by each P7 row and P5 column listed for
  each sample. For example, if sample_A has P7 rows 'A,D'
  and P5 columns '5,3', expand_sample_rows returns a list
  that includes two rows for sample_A where the first sample
  row has A for the P7 row and 5 for the P5 column and the
  second row has D and 3 for the P7 row and P5 column,
  respectively.
  """   
  # If P7 format != p7_rows or P5 format != p5_columns, then
  # return without expanding.
  pcr7_column, pcr5_column = test_pcr_format(column_name_list)
  if(pcr7_column == None or pcr5_column == None): 
    return(samplesheet_row_list)
  
  # Expand samples by pairing the i-th P7 row with the i-th P5
  # column for each of the specified rows and columns. Each sample
  # row has one PCR row and one PCR column.
  new_samplesheet_row_list = []
  num_element = len(column_name_list)
  for irow, samplesheet_row in enumerate(samplesheet_row_list):
    element_coordinates = [str( irow + 2 ), chr( pcr7_column + ord( 'A' ))]
    row_list = expand_rows(samplesheet_row[pcr7_column], element_coordinates = element_coordinates)
    element_coordinates = [str( irow + 2 ), chr( pcr5_column + ord( 'A' ))]
    column_list = expand_columns(samplesheet_row[pcr5_column], element_coordinates = element_coordinates)
    if(len(row_list) != len(column_list)):
      print('Error: number of rows is not equal to number of columns in sample sheet row %d' % (irow + 2))
      sys.exit(-1)
    num_expand = len(row_list)
    for iexpand in range(num_expand):
      element_list = []
      for ielem in range(num_element):
        if(ielem == pcr7_column):
          element_list.append(row_list[iexpand])
        elif(ielem == pcr5_column): 
          element_list.append(column_list[iexpand])
        else:
          element_list.append(samplesheet_row[ielem])
      new_samplesheet_row_list.append(element_list)
  return(new_samplesheet_row_list) 


def set_default_lanes_value( row_out_list, number_lanes ):
  #
  # Check that either all lanes are zero length strings
  # or they all are non-zero length.
  #
  zero_length_flag = False
  if( len( row_out_list[0]['lanes'] ) == 0 ):
    zero_length_flag = True
  for sample_dict in row_out_list:
    if( ( len( sample_dict['lanes'] ) >  0 and zero_length_flag == True) or
        ( len( sample_dict['lanes'] ) == 0 and zero_length_flag == False ) ):
      print('Error: mix of lanes with and without lane specifications', file=sys.stderr)
      sys.exit(-1)

  #
  # Set lanes to all lanes if zero_length_flag is True. Otherwise,
  # clean up lanes specs, if necessary.
  #
  if( zero_length_flag ):
    for sample_dict in row_out_list:
      sample_dict[ 'lanes' ] = list( range( 1, number_lanes + 1 ) )
  return( 0 )


# Finish by adding lanes and hash_file.
def make_samplesheet_indexes( column_name_list, samplesheet_row_list, number_lanes ):
  """
  Make well index lists for rt, p5, p7 barcode wells from the input samplesheet information.
  """
  num_col = len( column_name_list )
  row_out_list = []
  for irow, row_elements in enumerate( samplesheet_row_list ):
    if( len( row_elements ) < num_col ):
      print( 'Error: missing cells in row %d: %s' % ( irow + 1, ', '.join('"{0}"'.format(e) for e in row_elements ) ), file=sys.stderr )
      sys.exit( -1 )
    icol = 0
    # Initialize optional column values.
    external_sample_name = ''
    tissue = ''
    wrap_group = ''
    hash_file = ''
    sample_flags_list = ''
    library = ''
    process_group = ''
    rt_file = default_rt_file
    ligation_file = default_ligation_file
    p7_file = default_p7_file
    p5_file = default_p5_file
    lane_list = []
    for element_string, column_name_dict in zip( row_elements, column_name_list ):
      icol += 1

      element_coordinates = [ str( irow + 2 ), chr( icol + ord( 'A' ) - 1 ) ]
      if( column_name_dict['type'] == 'rt' ):
        rt_index_list = parse_wells( element_string, True, element_coordinates )
      elif( column_name_dict['type'] == 'p7' ):
        if( column_name_dict['format'] == 'wells' ):
          p7_index_list = parse_wells( element_string, True, element_coordinates )
        elif( column_name_dict['format'] == 'rows' ):
          p7_index_list = parse_rows( element_string, element_coordinates )
        else:
          print( 'Error: unexpected P7 format', file=sys.stderr )
          sys.exit( -1 )
      elif( column_name_dict['type'] == 'p5' ):
        if( column_name_dict['format'] == 'wells' ):
          p5_index_list = parse_wells( element_string, False, element_coordinates )
        elif( column_name_dict['format'] == 'columns' ):
          p5_index_list = parse_columns( element_string, element_coordinates )
        else:
          print( 'Error: unexpected P5 format', file=sys.stderr )
          sys.exit( -1 )
      elif( column_name_dict['type'] == 'sample_name' ):
          sample_name = element_string
      elif( column_name_dict['type'] == 'genome' ):
          genome = element_string
      elif( column_name_dict['type'] == 'external_sample_name' ):
          external_sample_name = element_string
      elif( column_name_dict['type'] == 'tissue' ):
          tissue = element_string
      elif( column_name_dict['type'] == 'wrap_group' ):
          wrap_group = element_string
      elif( column_name_dict['type'] == 'sample_flags' ):
        sample_flags_list = parse_sample_flags( element_string, element_coordinates)
      elif( column_name_dict['type'] == 'lanes' ):
        lane_list = parse_lanes( element_string, element_coordinates)
      elif( column_name_dict['type'] == 'rt_file' ):
          rt_file = element_string
      elif( column_name_dict['type'] == 'ligation_file' ):
          ligation_file = element_string
      elif( column_name_dict['type'] == 'p7_file' ):
          p7_file = element_string
      elif( column_name_dict['type'] == 'p5_file' ):
          p5_file = element_string
      elif( column_name_dict['type'] == 'library' ):
          library = element_string
      elif( column_name_dict['type'] == 'process_group' ):
          process_group = element_string
      elif( column_name_dict['type'] == 'hash_file' ):
          hash_file = element_string


    #
    row_out_list.append( { 'sample_name': sample_name,
                           'rt_index_list': rt_index_list,
                           'p7_index_list': p7_index_list,
                           'p5_index_list': p5_index_list,
                           'genome': genome,
                           'external_sample_name': external_sample_name,
                           'tissue': tissue,
                           'wrap_group': wrap_group,
                           'sample_flags': sample_flags_list,
                           'lanes' : lane_list,
                           'hash_file' : hash_file,
                           'rt_file' : rt_file,
                           'ligation_file' : ligation_file,
                           'p7_file' : p7_file,
                           'p5_file' : p5_file,
                           'library' : library,
                           'process_group' : process_group } )

  #
  # Set default lanes value if necessary.
  #
  set_default_lanes_value( row_out_list = row_out_list, number_lanes = number_lanes )

  return( row_out_list )


# OK
def check_pcr_indexes( row_out_list ):
  for irow, row_out in enumerate( row_out_list ):
    if( row_out['p7_index_list'] == [0] and row_out['p5_index_list'] == [0] ):
      print( 'Error: row %d: both p7 and p5 well indices are 0' % ( irow + 2 ) )
      sys.exit( -1 )

# OK
def dump_row_out_list( row_out_list ):
  """
  Diagnostic function to dump (barcode) well index lists.
  """
  print()
  print('** Report wells calculated from well indices to check indices as diagnostic. **')
  print()
  for row_out in row_out_list:
    print( 'sample_name: %s' % ( row_out['sample_name'] ) )
    print( '  rt_index_list: %s  wells: %s' % ( make_index_string( row_out['rt_index_list'] ), index_string_to_well_string( make_index_string( row_out['rt_index_list'] ), across_row_first=True, show_plate=True ) ) )
    print( '  p7_index_list: %s  wells: %s' % ( make_index_string( row_out['p7_index_list'] ), index_string_to_well_string( make_index_string( row_out['p7_index_list'] ), across_row_first=True, show_plate=False ) ) )
    print( '  p5_index_list: %s  wells: %s' % ( make_index_string( row_out['p5_index_list'] ), index_string_to_well_string( make_index_string( row_out['p5_index_list'] ), across_row_first=False, show_plate=False ) ) )
    print( '  lanes:  %s' % ( make_index_string( row_out['lanes']) ) )
    print( '  process group: %s' % ( row_out['process_group'] ) )
    print( '  genome: %s' % ( row_out['genome'] ) )
  return( 0 )


# OK
def get_pcr_row_col( column_name_list, samplesheet_row_list ):
  """
  Given lists of PCR rows and columns by sample, return lists of
  distinct PCR rows and columns for JSON output file. These values
  are used by the demux_dash.
  """ 
  # find required samplesheet column
  p5_samplesheet_col = None
  p7_samplesheet_col = None
  for icol, column_name_dict in enumerate( column_name_list ):
    if( column_name_dict['type'] == 'p5' ):            
      p5_samplesheet_col = icol                        
    if( column_name_dict['type'] == 'p7' ):
      p7_samplesheet_col = icol
                                
  # gather plate values         
  pair_samplesheet_row_list = []
  for samplesheet_row in samplesheet_row_list:
    p5_samplesheet_row_list = []
    for value_range in samplesheet_row[p5_samplesheet_col].split( ',' ):
      if(value_range == 'none'):
        p5_samplesheet_row_list.append( '0' )
        continue
      value_range = re.sub( r'\s', '', value_range )
      mobj = re.match( p5_re_col_pattern_2, value_range.strip() )
      if( not mobj ):
        print( 'Error: bad value or value range \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )
      col1 = int( mobj.group( 1 ) )
      col2 = col1
      if( mobj.group( 2 ) ):
        col2 = int( mobj.group( 3 ) )
      for col in range( col1, col2 + 1 ):
        p5_samplesheet_row_list.append( str( col ) )

    p7_samplesheet_row_list = []
    for value_range in samplesheet_row[p7_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      mobj = re.match( p7_re_row_pattern_2, value_range.strip() )
      if( not mobj ):
        print( 'Error: bad value or value range \'%s\'' % ( value_range ), file=sys.stderr )
        sys.exit( -1 )
      row1 = mobj.group( 1 )
      row2 = row1
      if( mobj.group( 2 ) ):
        row2 = mobj.group( 3 )
      for orow in range( ord( row1 ), ord( row2 ) + 1 ):
        p7_samplesheet_row_list.append( chr( orow ) )

    pair_samplesheet_row_list.append( [ p5_samplesheet_row_list, p7_samplesheet_row_list ] )

  pair_list_dict = {}
  for pair_samplesheet_row in pair_samplesheet_row_list:
    key = '_'.join( pair_samplesheet_row[0] ) + '_' + '_'.join( pair_samplesheet_row[1] )
    pair_list_dict.setdefault( key, pair_samplesheet_row )

  p5_col_list = []
  p7_row_list = []
  for key in pair_list_dict.keys():
    p5_col_list.extend( pair_list_dict[key][0] )
    p7_row_list.extend( pair_list_dict[key][1] )

  return( p5_col_list, p7_row_list )


# OK maybe
# bge
def get_pcr_wells( column_name_list, samplesheet_row_list ):
  """
  Given the samplesheet rows, return lists of specified PCR wells,
  if any.
  """
  p5_samplesheet_col = None
  p7_samplesheet_col = None
  for icol, column_name_dict in enumerate( column_name_list ):
    if( column_name_dict['type'] == 'p5' ):
      p5_samplesheet_col = icol
    if( column_name_dict['type'] == 'p7' ):
      p7_samplesheet_col = icol

  pair_samplesheet_row_list = []
  for samplesheet_row in samplesheet_row_list:
    p5_samplesheet_row_list = []
    for value_range in samplesheet_row[p5_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      if( mobj := re.match( r'([pP][0]?([1-4*])[-])?([a-hA-H])([0]?[1-9][0-2]?)([:]([pP][0]?([1-4*])[-])?([a-hA-H])([0]?[1-9][0-2]?))?$', value_range.strip() ) ):
        # Assume that the well-range strings were checked earlier in parse_wells() call so
        # those tests are omitted here.
        # first well
        row1 = mobj.group( 3 )
        col1 = int( mobj.group( 4 ) )
        row2 = row1
        col2 = col1
        if( col1 < 1 or col1 > 12 ):
          print( 'Error: bad well: \'%s\'' % ( value_range ), file=sys.stderr )
          sys.exit( -1 )
        if( mobj.group( 5 ) ):
          # second well, if this is a range
          row2 = mobj.group( 8 )
          col2 = int( mobj.group( 9 ) )
        across_row_first = False # p5
        index1 = well_to_index( 1, row1, col1, across_row_first, [ None, None ] )
        index2 = well_to_index( 1, row2, col2, across_row_first, [ None, None ] )
        for well_index in range( index1, index2 + 1 ):
          ipl, well = index_to_well( well_index - 1, across_row_first )
          p5_samplesheet_row_list.append( well )
      elif( mobj := re.match( r'0$', value_range.strip() ) ):
        p5_samplesheet_row_list.append( 'none' )
      else:
        if( not mobj ):
          print( 'Error: bad well or well range \'%s\'' % ( value_range ), file=sys.stderr )
          sys.exit( -1 )

    p7_samplesheet_row_list = []
    for value_range in samplesheet_row[p7_samplesheet_col].split( ',' ):
      value_range = re.sub( r'\s', '', value_range )
      if( mobj := re.match( r'([pP][0]?([1-4*])[-])?([a-hA-H])([0]?[1-9][0-2]?)([:]([pP][0]?([1-4*])[-])?([a-hA-H])([0]?[1-9][0-2]?))?$', value_range.strip() ) ):
        # Assume that the well-range strings were checked earlier in parse_wells() call so
        # those tests are omitted here.
        # first well
        row1 = mobj.group( 3 )
        col1 = int( mobj.group( 4 ) )
        row2 = row1
        col2 = col1
        if( col1 < 1 or col1 > 12 ):
          print( 'Error: bad well: \'%s\'' % ( value_range ), file=sys.stderr )
          sys.exit( -1 )
        if( mobj.group( 5 ) ):
          # second well, if this is a range
          row2 = mobj.group( 8 )
          col2 = int( mobj.group( 9 ) )
        across_row_first = True # p7
        index1 = well_to_index( 1, row1, col1, across_row_first, [ None, None ] )
        index2 = well_to_index( 1, row2, col2, across_row_first, [ None, None ] )
        for well_index in range( index1, index2 + 1 ):
          ipl, well = index_to_well( well_index - 1, across_row_first )
          p7_samplesheet_row_list.append( well )
      elif( mobj := re.match( r'0$', value_range.strip() ) ):
        p7_samplesheet_row_list.append( 'none' )
      else:
        if( not mobj ):
          print( 'Error: bad well or well range \'%s\'' % ( value_range ), file=sys.stderr )
          sys.exit( -1 )
    pair_samplesheet_row_list.append( [ p5_samplesheet_row_list, p7_samplesheet_row_list ] )

  pair_list_dict = {}
  for pair_samplesheet_row in pair_samplesheet_row_list:
    key = '_'.join( pair_samplesheet_row[0] ) + '_' + '_'.join( pair_samplesheet_row[1] )
    pair_list_dict.setdefault( key, pair_samplesheet_row )

  p5_well_list = []
  p7_well_list = []
  for key in pair_list_dict.keys():
    p5_well_list.extend( pair_list_dict[key][0] )
    p7_well_list.extend( pair_list_dict[key][1] )

  return( p5_well_list, p7_well_list )


def write_samplesheet_json_format( file, column_name_list, samplesheet_row_list, row_out_list, wrap_groups_dict, level = 3, number_lanes = 1, sequencer_run_directory = 'NA', sequencer_class = 'illumina' ):
  """
  Write an output samplesheet file in JSON format.
  """
  # Store input samplesheet for for reference if questions arise.
  input_samplesheet_rows = []
  for samplesheet_row in samplesheet_row_list:
    input_samplesheet_rows.append( ','.join( '"{0}"'.format( e ) for e in samplesheet_row ) )

  # Store sample information for processing pipeline.
  sample_index_list = []
  for row_out in row_out_list:
    if( row_out['lanes'] != '' ):
      lanes = make_index_string( row_out['lanes'] )
    else:
      lanes = ''
    if( row_out['process_group'] != '' ):
      process_group = row_out['process_group']
    else:
      process_group = '1'

    sample_index_list.append( { 'sample_id' : row_out['sample_name'],
                                'ranges' : ':'.join( [ make_index_string( row_out['rt_index_list'] ),
                                                       make_index_string( row_out['p7_index_list'] ),
                                                       make_index_string( row_out['p5_index_list'] ) ] ),
                                'lanes' : lanes,
                                'tissue' : row_out['tissue'],
                                'genome' : row_out['genome'],
                                'hash_file' : row_out['hash_file'],
                                'sample_flags' : ','.join( row_out['sample_flags'] ),
                                'external_sample_name' : row_out['external_sample_name'],
                                'wrap_group' : row_out['wrap_group'],
                                'rt_file' : row_out['rt_file'],
                                'ligation_file' : row_out['ligation_file'],
                                'p7_file' : row_out['p7_file'],
                                'p5_file' : row_out['p5_file'],
                                'library' : row_out['library'],
                                'process_group' : process_group })

  # Store information for dashboard(s).

  # Note: assume that the header was checked for consistent
  #       p5 => columns and p7 => rows
  pcr_format = None
  for icol, column_name_dict in enumerate( column_name_list ):
    if( column_name_dict['type'] == 'p5' ):
      if( column_name_dict['format'] == 'columns' ):
        pcr_format = 'row_col'
      elif( column_name_dict['format'] == 'indexes' ):
        pcr_format = 'indexes'
      elif( column_name_dict['format'] == 'wells' ):
        pcr_format = 'wells'

  # PCR rows and columns specified?
  p5_col_list = None
  p7_row_list = None
  if( pcr_format == 'row_col' ):
    p5_col_list, p7_row_list = get_pcr_row_col( column_name_list, samplesheet_row_list )

  # PCR wells specified?
  p5_well_list = None
  p7_well_list = None
  if( pcr_format == 'wells' ):
    p5_well_list, p7_well_list = get_pcr_wells( column_name_list, samplesheet_row_list )

  # Tissue dict.
  tissue_dict = {}
  for row_out in row_out_list:
    tissue_dict[row_out['sample_name']] = row_out.get('tissue', '')

  # external_sample_name dict.
  external_sample_name_dict = {}
  for row_out in row_out_list:
    external_sample_name_dict[row_out['sample_name']] = row_out.get('external_sample_name', '')

#  # Make input column name list.
#  input_column_name_list = []
#  for column_name_dict in column_name_list:
#    if(column_name_dict.get('format', None) == None):
#      input_column_name_list.append(column_name_dict['type'])
#    else:
#      input_column_name_list.append(column_name_dict['type'] + '_' + column_name_dict['format'])
#  input_column_names = ','.join(input_column_name_list)

  # Make input column name list.
  input_column_name_list = []
  for column_name_dict in column_name_list:
    input_column_name_list.append(column_name_dict)

  # JSON structure.
  sample_data = { 'json_file_version' : json_file_version,
                  'sequencer_class' : sequencer_class,
                  'sequencer_run_directory' : sequencer_run_directory,
                  'level' : level,
                  'number_lanes' : number_lanes,
                  'input_samplesheet_column_names' : input_column_name_list,
                  'input_samplesheet_rows' : input_samplesheet_rows,
                  'p5_col_list' : p5_col_list,
                  'p7_row_list' : p7_row_list,
                  'p5_well_list' : p5_well_list,
                  'p7_well_list' : p7_well_list,
                  'sample_index_list' : sample_index_list
                }

  file.write(json.dumps(sample_data, indent=4))

  return( 0 )


#
# Count distinct well indices.
#
def count_wells( index_list ):
  num_well = len( set( index_list ) )
  return( num_well )


#
# Gather wrap group information. That is, the
# wrap group names, if they are given in the
# input samplesheet file, and the samples that belong
# to them.
#
def get_wrap_groups( row_out_list ):
  wrap_groups_dict = {}
  for irow, row_out in enumerate(row_out_list):
    # Gather wrap_group values in wrap_groups_dict.
    # Does this sample have a wrap group value?
    if( row_out.get( 'wrap_group' ) != None and row_out['wrap_group'] != '' ):
      if( wrap_groups_dict.get( row_out['wrap_group'] ) == None ):
        wrap_groups_dict[row_out['wrap_group']] = [ row_out['sample_name'] ]
      elif( row_out['sample_name'] not in wrap_groups_dict[row_out['wrap_group']] ):
        wrap_groups_dict[row_out['wrap_group']].append( row_out['sample_name'] )
  return( wrap_groups_dict )
  
  
    
def samplesheet_report( samplesheet_row_list, row_out_list, wrap_groups_dict, args ):
  print()
  print( '== Samplesheet information ==' )
  print( '  Level:             %s' % ( args.level ) )
 
  # Report sample names after replacing unacceptable characters.
  string_list = []
  print( '  Sample names after converting unacceptable characters to \'.\':' )
  for irow, row_out in enumerate(row_out_list):
    if( row_out_list[irow]['sample_name'] in string_list):
      continue
    print( '    %s' % ( row_out['sample_name'] ) )
    string_list.append(row_out_list[irow]['sample_name'])

  # Sort rows by lane and sample_name cells.
  row_out_list = sorted(row_out_list, key = lambda x: (x['lanes'], x['sample_name']))

  # Calculate spacing.
  max_len_samplename = 0
  for row_out in row_out_list:
    if( len( row_out['sample_name'] ) > max_len_samplename ):
      max_len_samplename = len( row_out['sample_name'] )
  if( len( 'name' ) > max_len_samplename ):
    max_len_samplename = len( 'name' )
  nspace_samplename = max_len_samplename - len( 'name' )

  lane_set_dict = dict()
  max_len_laneset = 0
  for row_out in row_out_list:
    lane_set_dict.setdefault(make_index_string(row_out['lanes']), 0)
    if(len(make_index_string(row_out['lanes'])) > max_len_laneset):
      max_len_laneset = len(make_index_string(row_out['lanes']))
  if( len( 'lanes' ) > max_len_laneset ):
    max_len_laneset = len( 'lanes' )
  nspace_laneset = max_len_laneset - len( 'lanes' )

  # Count lane sets.
  num_laneset = len(lane_set_dict.keys())

  # Report well counts by lane and sample.
  print( '  Sample well counts:' )
  if( num_laneset == 1 ):
    print( '    name%s    RT    P7    P5' % ( ' ' *  nspace_samplename ) )
  else:
    print( '    lanes%s    name%s    RT    P7    P5' % ( ' ' * nspace_laneset, ' ' * nspace_samplename ) )

  for irow, row_out in enumerate(row_out_list):
    if(irow > 0
       and row_out_list[irow]['lanes'] == row_out_list[irow-1]['lanes']
       and row_out_list[irow]['sample_name'] == row_out_list[irow-1]['sample_name']
       and count_wells(row_out_list[irow]['rt_index_list']) == count_wells(row_out_list[irow-1]['rt_index_list'])
       and count_wells(row_out_list[irow]['p7_index_list']) == count_wells(row_out_list[irow-1]['p7_index_list'])
       and count_wells(row_out_list[irow]['p5_index_list']) == count_wells(row_out_list[irow-1]['p5_index_list'])):
      continue

    if( num_laneset == 1 ):
      print( '    %s%s    %d    %d    %d' % ( row_out['sample_name'],
                                            ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                            count_wells( row_out['rt_index_list'] ),
                                            count_wells( row_out['p7_index_list'] ),
                                            count_wells( row_out['p5_index_list'] ) ) )
    else:
      print( '    %s%s    %s%s    %d    %d    %d' % ( make_index_string(row_out['lanes']),
                                                    ' ' * (  max_len_laneset - len( make_index_string( row_out['lanes'] ) ) ),
                                                    row_out['sample_name'],
                                                    ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                                    count_wells( row_out['rt_index_list'] ),
                                                    count_wells( row_out['p7_index_list'] ),
                                                    count_wells( row_out['p5_index_list'] ) ) )

  # Report sample genome names.
  print( '  Sample genomes:' )
  max_len_genome = 0
  for row_out in row_out_list:
    if( len( row_out['genome'] ) > max_len_genome ):
      max_len_genome = len( row_out['genome'] )
  if( len( 'genome' ) > max_len_genome ):
    max_len_genome = len( 'genome' )
  print( '    name%s    genome' % ( ' ' * ( max_len_samplename - len( 'name' ) ) ) )

  for irow, row_out in enumerate(row_out_list):
    if(irow > 0
       and row_out_list[irow]['sample_name'] == row_out_list[irow-1]['sample_name']
       and row_out_list[irow]['genome'] == row_out_list[irow-1]['genome']):
      continue

    print( '    %s%s    %s' % ( row_out['sample_name'],
                                          ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                          row_out['genome'] ) )

  # Report sample flags
  sample_flags_flag = False
  for row_out in row_out_list:
    if( 'sample_flags' in row_out and len( ''.join( row_out['sample_flags'] ) ) > 0 ):
      sample_flags_flag = True
  if( sample_flags_flag == True ):
    max_len_samplename = 0
    for row_out in row_out_list: 
      if( len( row_out['sample_name'] ) > max_len_samplename ):
        max_len_samplename = len( row_out['sample_name'] )
    print( '  Sample flags:' )

    if( num_laneset == 1 ):
      print( '    name%s    %s' % ( ' ' * nspace_samplename, 'sample_flags' ) )
    else:
      print( '    lanes%s    name%s    %s' % ( ' ' * nspace_laneset, ' ' * nspace_samplename, 'sample_flags' ) )

    for irow, row_out in enumerate( row_out_list ):
      if( irow == 0 or
          make_index_string( row_out_list[irow]['lanes'] ) != make_index_string( row_out_list[irow-1]['lanes'] ) or
          row_out_list[irow]['sample_flags'] != row_out_list[irow-1]['sample_flags'] or
          ''.join( row_out_list[irow]['sample_flags'] ) != ''.join( row_out_list[irow-1]['sample_flags'] ) ):
        if( ''.join(row_out['sample_flags']) != '' ):
          if( num_laneset == 1 ):
            print('    %s%s    %s' % ( row_out['sample_name'], ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ), ', '.join( row_out['sample_flags'] ) ) )
          else:
            print('    %s%s    %s%s    %s' % ( make_index_string( row_out['lanes'] ), ' ' * ( max_len_laneset - len( make_index_string( row_out['lanes'] ) ) ),
                                               row_out['sample_name'],  ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                               ', '.join(row_out['sample_flags'] ) ) )

  # Report external sample names.
  external_sample_name_flag = False
  for row_out in row_out_list:
    if( 'external_sample_name' in row_out and len( row_out['external_sample_name'] ) > 0 ):
      external_sample_name_flag = True
  if( external_sample_name_flag == True ):
    print( '  External sample names with tabs converted to spaces:' )
    if( num_laneset == 1 ):
      print( '    name%s    external_sample_name' % ( ' ' *  nspace_samplename ) )
    else:
      print( '    lanes%s    name%s    external_sample_name' % ( ' ' * nspace_laneset, ' ' * nspace_samplename ) )

    for irow, row_out in enumerate( row_out_list ):
      if( irow == 0 or
          make_index_string( row_out_list[irow]['lanes'] ) != make_index_string( row_out_list[irow-1]['lanes'] ) or
          row_out_list[irow]['sample_name'] != row_out_list[irow-1]['sample_name'] or
          row_out_list[irow]['external_sample_name'] != row_out_list[irow-1]['external_sample_name'] ):
        if( num_laneset == 1 ):
          print( '    %s%s     %s' % ( row_out['sample_name'],
                                       ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                       row_out.get('external_sample_name', '' ) ) )
        else:
          print( '    %s%s    %s%s    %s' % ( make_index_string( row_out['lanes'] ),
                                              ' ' * ( max_len_laneset - len( make_index_string( row_out['lanes'] ) ) ),
                                              row_out['sample_name'],
                                              ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                              row_out.get('external_sample_name', '' ) ) )
 
  # Report tissue names.
  tissue_flag = False
  for row_out in row_out_list:
    if( 'tissue' in row_out and len( row_out['tissue'] ) > 0 ):
      tissue_flag = True
  if( tissue_flag == True ):
    print( '  Tissue assignments with tabs converted to spaces:' )
    if( num_laneset == 1 ):
      print( '    name%s    tissue' % ( ' ' *  nspace_samplename ) )
    else:
      print( '    lanes%s    name%s    tissue' % ( ' ' * nspace_laneset, ' ' * nspace_samplename ) )
    for irow, row_out in enumerate(row_out_list):
      if( irow == 0 or
          make_index_string( row_out_list[irow]['lanes'] ) != make_index_string( row_out_list[irow-1]['lanes'] ) or
          row_out['sample_name'] != row_out_list[irow-1]['sample_name'] or
          row_out['tissue'] != row_out_list[irow-1]['tissue'] ):
        if( num_laneset == 1 ):
          print('    %s%s    %s' % ( row_out['sample_name'],
                                     ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                     row_out.get('tissue', '')))
        else:
          print('    %s%s    %s%s    %s' % ( make_index_string( row_out['lanes'] ),
                                              ' ' * ( max_len_laneset - len( make_index_string( row_out['lanes'] ) ) ),
                                             row_out['sample_name'],
                                             ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                             row_out.get('tissue', '')))

  # Report information about wrap groups for the wrapping,
  # if it exists in the input samplesheet file.
  wrap_group_flag = False
  for row_out in row_out_list:
    if( 'wrap_group' in row_out and len( row_out['wrap_group'] ) > 0 ):
      wrap_group_flag = True
  if( wrap_group_flag == True ):
    print( '  Wrap group assignments:' )
    if( num_laneset == 1 ):
      print( '    name%s    wrap_group' % ( ' ' *  nspace_samplename ) )
    else:
      print( '    lanes%s    name%s    wrap_group' % ( ' ' * nspace_laneset, ' ' * nspace_samplename ) )
    for irow, row_out in enumerate(row_out_list):
      if( irow == 0 or
          make_index_string( row_out_list[irow]['lanes'] ) != make_index_string( row_out_list[irow-1]['lanes'] ) or
          row_out['sample_name'] != row_out_list[irow-1]['sample_name'] or
          row_out['wrap_group'] != row_out_list[irow-1]['wrap_group'] ):
        if( num_laneset == 1 ):
          print('    %s%s    %s' % ( row_out['sample_name'],
                                     ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                     row_out.get('wrap_group', '')))
        else:
          print('    %s%s    %s%s    %s' % ( make_index_string( row_out['lanes'] ),
                                              ' ' * ( max_len_laneset - len( make_index_string( row_out['lanes'] ) ) ),
                                             row_out['sample_name'],
                                             ' ' * ( max_len_samplename - len( row_out['sample_name'] ) ),
                                             row_out.get('wrap_group', '')))

  print( '  Illumina run directory: %s' % ( args.run_dir ) )
  print( '  Run scirna_samplesheet.py -d for more information.' )
  return( 0 )


def write_samplesheet_template():
  filename = 'samplesheet.template.csv'
  with open( filename, 'wt' ) as fp:
    print( 'rt_wells,p7_rows,p5_columns,sample_name,genome,external_sample_name,tissue,wrap_group,lanes,hash_file,sample_flags', file=fp )
  return( 0 )


if __name__ == '__main__':
  parser = argparse.ArgumentParser(description='A program to convert sci-RNA CSV samplesheet to pipeline samplesheet.')
  parser.add_argument('-i', '--input', required=False, default=None, help='Input CSV samplesheet filename (required string).')
  parser.add_argument('-o', '--output', required=False, default=None, help='Output samplesheet filename (required string).')
  parser.add_argument('-s', '--sequencer_class', required=False, default='illumina', help='Sequencing machine class (not required string: \'illumina\' or \'ultima\'. Default: illumina')
  parser.add_argument('-r', '--run_dir', required=False, default=None, help='Illumina run directory path (optional string).')
  parser.add_argument('-l', '--level', type=int, required=False, choices=[ 2, 3 ], default=3, help='Two or three level sci-RNA-seq experiment (default: %(default)s) (optional integer).')
  parser.add_argument('-n', '--number_lanes', type=int, required=False, default=None, help='Number of flowcell lanes (required integer).')
  parser.add_argument('-e', '--template', required=False, action='store_true', help='Write template samplesheet file (\'samplesheet.template.csv\') with standard column formats and exit (optional flag).')
  parser.add_argument('-d', '--documentation', required=False, action='store_true', help='Display documentation and exit (optional flag).')
  parser.add_argument('-v', '--version', action='version', version=program_version)
  args = parser.parse_args()

  #
  # Need to do messages.
  #
  print('================================================================================', file=sys.stderr)
  print('Need to do:', file=sys.stderr)
  print('  o  add a command line option to dump diagnostic information', file=sys.stderr)
  print('================================================================================', file=sys.stderr)
  print('', file=sys.stderr)


  # Write documentation.
  if( args.documentation ):
    display_documentation()
    sys.exit( 0 )

  # Write samplesheet template file.
  if( args.template ):
    write_samplesheet_template()
    sys.exit( 0 )

  # Check for required command line parameters.
  error_string = '' 
  if( args.input == None ):
    error_string += '  input filename parameter: -i <input filename> or --input <input filename>\n'
  if( args.output == None ):
    error_string += '  output filename parameter: -o <output filename> or --output <output filename>\n'
  if( args.number_lanes == None):
    error_string += '  number_lanes parameter: -n <number_lanes>\n'
  if( len( error_string ) > 0 ):
    print( 'Error: missing command line parameters\n%s' % ( error_string ) )
    print( 'For help run \'scirna_samplesheet.py -h\' or \'scirna_samplesheet.py -d\'' )
    sys.exit( -1 )

  #
  # Check command line parameter consistency.
  #
  check_args( args )
  
  # Go to work.
  filename_in = args.input
  filename_out = args.output
  
  column_name_list, samplesheet_row_list = read_samplesheet( open( filename_in, newline='' ) )

  clean_samplesheet_data( column_name_list, samplesheet_row_list )
  print('== Check input samplesheet.')
  samplesheet_row_list = check_sample_names( column_name_list, samplesheet_row_list )
  check_genome_names( column_name_list, samplesheet_row_list )
  check_external_sample_name( column_name_list, samplesheet_row_list )
  check_tissue( column_name_list, samplesheet_row_list )
  check_wrap_group( column_name_list, samplesheet_row_list )
  check_lanes( column_name_list, samplesheet_row_list )
  check_hash_file( column_name_list, samplesheet_row_list )
  check_barcode_files( column_name_list, samplesheet_row_list )
  check_lane_sample_consistency( column_name_list, samplesheet_row_list )
  check_process_groups(column_name_list, samplesheet_row_list )

#  if(not args.no_expand_pcr_rows_columns):
  samplesheet_row_list = expand_sample_rows(column_name_list, samplesheet_row_list)

  row_out_list = make_samplesheet_indexes( column_name_list = column_name_list, samplesheet_row_list = samplesheet_row_list, number_lanes = args.number_lanes )
  check_pcr_indexes( row_out_list )
  wrap_groups_dict = get_wrap_groups( row_out_list )
  write_samplesheet_json_format( open( filename_out, 'w' ), column_name_list, samplesheet_row_list, row_out_list, wrap_groups_dict, level = args.level, number_lanes = args.number_lanes, sequencer_run_directory = args.run_dir, sequencer_class = args.sequencer_class )
  samplesheet_report( samplesheet_row_list, row_out_list, wrap_groups_dict, args )
  # diagnostic dump
  dump_row_out_list( row_out_list )


