#!/bin/bash

NOW=`date '+%Y%m%d_%H%M%S'`
WORK_DIR="$PWD/work_demux"
TRACE_FILE="$PWD/trace.demux.${NOW}.tsv"
CONFIG_FILE="$PWD/experiment.config"

/net/gs/vol1/home/bge/bin/nextflow run /net/gs/vol1/home/bge/git/bbi-scirna-demux/main.nf -c $CONFIG_FILE -w $WORK_DIR -with-trace $TRACE_FILE  -resume
