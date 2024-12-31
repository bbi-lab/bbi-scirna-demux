#!/bin/bash

ns=$1

../scirna_samplesheet.py -i samplesheet${ns}.csv -o foo.json -n 4
echo
cat samplesheet${ns}.csv
