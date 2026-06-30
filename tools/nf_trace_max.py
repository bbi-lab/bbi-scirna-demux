#!/usr/bin/env python3

import sys
import re
import csv


#       0    1            2          3         4     5         6         7               8           9          10               11          12      1 
# task_id hash    native_id       name    status  exit    submit  duration        realtime        %cpu    peak_rss        peak_vmem       rchar   wchar
# 2       4f/e74652       8255757 make_sample_map_json (1)        COMPLETED       0       2026-02-26 14:52:32.174 4.2s    284ms   41.8%   11 MB   25.7 MB 894.9 KB        3.3 KB
# 1       11/0ad851       8255756 make_merge_demux_json (1)       COMPLETED       0       2026-02-26 14:52:32.096 4.3s    153ms   77.6%   3.1 MB  4.3 MB  895.9 KB        32.5 KB
# 11      90/e8e4ed       8255764 merge_demux (9) COMPLETED       0       2026-02-26 14:52:36.970 4.5s    69ms    57.7%   3.1 MB  4.4 MB  177.8 KB        11.6 KB
# 13      07/f000a5       8255765 merge_demux (11)        COMPLETED       0       2026-02-26 14:52:37.046 4.5s    53ms    84.2%   3.1 MB  4.4 MB  167.7 KB        6.5 KB
# ...
# 5587    a5/4bebfe       6244286 align_bams (1527)       COMPLETED       0       2026-01-05 23:14:25.731 18m 25s 8m 35s  471.0%  14 GB   14.5 GB 27.3 GB 10 GB
# 5595    c5/1083a5       6244294 align_bams (1535)       COMPLETED       0       2026-01-05 23:15:45.706 17m 55s 7m 50s  452.6%  14 GB   14.5 GB 28.4 GB 10.7 GB
# 5591    33/763b9c       6244287 align_bams (1531)       COMPLETED       0       2026-01-05 23:14:30.612 21m 25s 11m 14s 533.4%  14 GB   14.5 GB 30.5 GB 12.2 GB
# ...
# 211     a8/c4e5b0       616912  make_cds_raw (1)        COMPLETED       0       2025-10-23 18:29:46.634 2h 29m 55s      2h 29m 48s      77.2%   13.4 GB 14.5 GB 2.3 TB  47.6 GB
# 214     3d/7d7367       618408  make_barnyard_json (1)  COMPLETED       0       2025-10-23 20:59:42.263 4.6s    275ms   36.9%   11 MB   25.7 MB 882.7 KB        1.6 KB
# 213     c8/9ec368       618409  assign_hash_raw (1)     FAILED  1       2025-10-23 20:59:42.384 4.5s    18ms    -       -       -       -       -
# 212     38/211ee4       618407  make_generate_qc_no_hash (1)    FAILED  0       2025-10-23 20:59:42.182 29.7s   211ms   21.0%   1.5 MB  4.3 MB  122.8 KB        1.2 KB
#


pobj_list_1 = list()
for i in range(5):
  new_list = [None, None]
  pobj_list_1.append(new_list)
pobj_list_1[0][0] = 0.001
pobj_list_1[0][1] = re.compile(r'([0-9.]+)ms')
pobj_list_1[1][0] = 1.0
pobj_list_1[1][1] = re.compile(r'([0-9.]+)s')
pobj_list_1[2][0] = 60.0
pobj_list_1[2][1] = re.compile(r'([0-9.]+)m')
pobj_list_1[3][0] = 3600.0
pobj_list_1[3][1] = re.compile(r'([0-9.]+)h')
pobj_list_1[4][0] = 86400.0
pobj_list_1[4][1] = re.compile(r'([0-9.]+)d')

def get_duration_in_seconds(string):
  parts = string.strip().split()
  # duration in seconds
  duration = 0.0
  for part in parts:
    for pobj in pobj_list_1:
      mobj = pobj[1].match(part)
      if(mobj != None):
        duration += float(mobj.group(1)) * pobj[0]
        break
  return(duration)


pobj_list_2 = list()
for i in range(4):
  new_list = [None, None]
  pobj_list_2.append(new_list)
pobj_list_2[0][0] = 1.0
pobj_list_2[0][1] = re.compile(r'([0-9.]+) B$')
pobj_list_2[1][0] = 1000.0
pobj_list_2[1][1] = re.compile(r'([0-9.]+) KB$')
pobj_list_2[2][0] = 1000000.0
pobj_list_2[2][1] = re.compile(r'([0-9.]+) MB$')
pobj_list_2[3][0] = 1000000000.0
pobj_list_2[3][1] = re.compile(r'([0-9.]+) GB$')

def get_mem_peak(string):
  mem_use = 0.0
  for pobj in pobj_list_2:
    mobj = pobj[1].match(string.strip())
    if(mobj != None):
      mem_use = float(mobj.group(1)) * pobj[0]
      break
  return(mem_use)


filename = sys.argv[1]
with open(filename, 'r', newline='') as ifh:
  csv_reader = csv.reader(ifh, delimiter='\t')
  process_dict = dict()
  header = next(csv_reader)
  for row in csv_reader:
    # Process name.
    process_name = row[3].split()[0]
    if(process_dict.get(process_name) == None):
      new_dict = dict()
      process_dict[process_name] = new_dict
    # Get maximum run duration
    if(process_dict[process_name].get('duration') == None):
      process_dict[process_name]['duration'] = 0.0
    duration_seconds = get_duration_in_seconds(row[7])
    if(duration_seconds > process_dict[process_name]['duration']):
      process_dict[process_name]['duration'] = duration_seconds
    # Get maximum vmem peak
    if(process_dict[process_name].get('vmem_peak') == None):
      process_dict[process_name]['vmem_peak'] = 0.0
    vmem_peak = get_mem_peak(row[11])
    if(vmem_peak > process_dict[process_name]['vmem_peak']):
      process_dict[process_name]['vmem_peak'] = vmem_peak
    # Get maximum rss_mem_peak
    if(process_dict[process_name].get('rss_mem_peak') == None):
      process_dict[process_name]['rss_mem_peak'] = 0.0
    rss_mem_peak = get_mem_peak(row[10])
    if(rss_mem_peak > process_dict[process_name]['rss_mem_peak']):
      process_dict[process_name]['rss_mem_peak'] = rss_mem_peak

  process_name_length_max = 0
  for process_name in process_dict.keys():
    if(len(process_name) > process_name_length_max):
      process_name_length_max = len(process_name)

  print('%-*s  %s  %s  %s' % (process_name_length_max, 'process name', 'max. duration (min.)', 'max vmem peak (MB)', 'max rss mem peak (MB)'))
  for process_name in process_dict.keys():
    print('%-*s   %12.3f          %12.1f           %12.1f' % (process_name_length_max, process_name, process_dict[process_name]['duration'] / 60.0, process_dict[process_name]['vmem_peak'] / 1000000.0, process_dict[process_name]['rss_mem_peak'] / 1000000.0))


