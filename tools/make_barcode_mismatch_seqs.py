#!/usr/bin/env python3

import sys


barcode_filename = '/net/gs/vol1/home/bge/git/bbi-scirna-demux/data/foo.out'

with open(barcode_filename) as fp:
  barcode_list = list()
  for barcode in fp:
    barcode_list.append(barcode.strip())

# print(barcode_list)
# print(len(barcode_list))

# Make mismatch sequences.
subs = list(['A', 'C', 'G', 'T', 'N'])
# subs = list(['A'])

for bc_wl in barcode_list:
  # Print original whitelist barcode.
#  print('bc_wl: %s' % (bc_wl))
  print('%s' % (bc_wl))
  # Generate single mismatch sequences including possible collisions.
  len_seq = len(bc_wl)
  for ibas in range(len_seq):
#    print('ibas: %d  %s' % (ibas+1, bc_wl[ibas:ibas+1]))
    for sbas in subs:
      if(sbas == bc_wl[ibas:ibas+1]):
        continue
      bc_mm = bc_wl[:ibas] + sbas + bc_wl[ibas+1:]
#      print(bc_wl)
      print(bc_mm)
#      print()
      

