#!/usr/bin/env python3

import argparse



if __name__ == '__main__':
  parser = argparse.ArgumentParser('Run scrublet.')
  parser.add_argument('-i', '--input', required=True, help='Input fastq filename.')
  args = parser.parse_args()

  in_filename=args.input
#   out_filename=''


# @25.0464-P5none-P7B01_331|25.0464|none|B01|P06-A02_LIG27|AGCTTAAG
# GGGTACATATCTCGGTCAGTGGCGGGATCCTTTGATAATGAAGGCATTGCCATTTTTGCGCTTCAGTTCACTTACTACTTATGG
# +
# IIIII-9I9IIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII9IIIIIIIIIIIIIIIIIIIIIIIIIIIII9III

#   ofh = open(out_filename, 'w'')

  with open(in_filename, 'r') as ifh:
  
    iline = 0
    for line in ifh:
      if(iline % 4 == 0):
        parts1 = line.strip().split('|')
        parts2 = parts1[0].split('_')
        header = '|'.join([parts2[0]]+parts1[1:])
      elif(iline % 4 == 1):
        print('%s %s' % (header, line.strip()))
#      if(iline > 20):
#        break 
      iline += 1
#   close(ofh)
