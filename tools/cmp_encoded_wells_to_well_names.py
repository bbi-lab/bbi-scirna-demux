#!/usr/bin/env python3

import sys
import string
import re


# 25.0464-P5none-P7P01-B01_1|25.0464|none|P01-B01|P06-A02_LIG58|AGCAGGGG

# CB:Z:AACTGAGAAAATGGAAAAATCAAAAAAA
# tag_bb = read.get_tag('bb')


def from_base(digits, base):
    """Converts a list of digits in the given base to an integer.

    The first digit is the most significant and the base is assumed to
    be an integer greater than or equal to 2.
    """
    power = 1
    number = 0
    for digit in reversed(digits):
        number += power * int(digit)
        power *= base
    return number


pobj1 = re.compile(r'P([01][0-9])-([A-H])([01][0-9])')
pobj2 = re.compile(r'none')

row_dict = {'A': 1, 'B': 2, 'C': 3, 'D': 4, 'E': 5, 'F': 6, 'G': 7, 'H': 8}
def well_to_index(well, along_row_first):
  if(pobj2.search(well) != None):
    return(0)

  mobj = pobj1.search(well)
  if(mobj == None):
    print('Error: unable to match well %s' % (well), file=sys.stderr)
    sys.exit(-1)
  plate = int(mobj.group(1))
  row = mobj.group(2)
  col = int(mobj.group(3))

  index = 0
  if(along_row_first == True):
    index = (plate-1) * 96
    index += (row_dict[row]-1) * 12
    index += col
  else:
   index = (plate-1) * 96
   index += (col-1) * 8
   index += row_dict[row]
  return(index)


in_filename = sys.argv[1]

with open(in_filename) as ifh:
  pobj3 = re.compile('LIG([0-9]+)')
  iline = 0
  btrans_tab = str.maketrans('ACGT', '0123')
  for line in ifh:
    iline += 1
    parts1 = line.strip().split()
    barcode = parts1[0]

    rteb4 = barcode[0:7]
    lgeb4 = barcode[7:14]
    p7eb4 = barcode[14:21]
    p5eb4 = barcode[21:]

    rtb4 = rteb4.translate(btrans_tab)
    lgb4 = lgeb4.translate(btrans_tab)
    p7b4 = p7eb4.translate(btrans_tab)
    p5b4 = p5eb4.translate(btrans_tab)
 
    rtidx = from_base(list(rtb4), 4)
    lgidx = from_base(list(lgb4), 4)
    p7idx = from_base(list(p7b4), 4)
    p5idx = from_base(list(p5b4), 4)

#    parts1 = query_name.strip().split('|')
#    parts2 = parts1[4].split('_')

    parts2 = parts1[1].split('_')
    p5swell = parts2[0]
    p7swell = parts2[1]
    rtswell = parts2[2]
    lgswell = parts2[3]

    rt_wi = well_to_index(rtswell, True)
    p7_wi = well_to_index(p7swell, True)
    p5_wi = well_to_index(p5swell, False)
    mobj = pobj3.search(lgswell)
    if(mobj == None):
      print('Error: unable to parse ligation barcode name')
      sys.exit(-1)
    lg_wi = int(mobj.group(1))

    if(rt_wi != rtidx):
      print('rt index mismatch', file=sys.stderr)
    if(lg_wi != lgidx):
      print('lg index mismatch', file=sys.stderr)
    if(p7_wi != p7idx):
      print('p7 index mismatch', file=sys.stderr)
    if(p5_wi != p5idx):
      print('p5 index mismatch', file=sys.stderr)


#    print('%s %s %s %s  %s %s %s %s  %s %s %s %s  %d %d %d %d  %d %d %d %d' % (parts2[0], parts2[1], parts1[3], parts1[2], rteb4, lgeb4, p7eb4, p5eb4, rtb4, lgb4, p7b4, p5b4, rtidx, lgidx, p7idx, p5idx, rt_wi, lg_wi, p7_wi, p5_wi))

    print('%s %s %s %s  %s %s %s %s  %s %s %s %s  %d %d %d %d  %d %d %d %d' % (rtswell, lgswell, p7swell, p5swell, rteb4, lgeb4, p7eb4, p5eb4, rtb4, lgb4, p7b4, p5b4, rtidx, lgidx, p7idx, p5idx, rt_wi, lg_wi, p7_wi, p5_wi))

#    if(iline > 10):
#      break



