process run_bclconvert {

  input:
    val samplesheet_json
    val illumina_run
    val p7_barcode_file
    val p5_barcode_file

  output:
    path "fastqs_bclconvert/*.fastq.gz"

  script:
  """
  module load modules modules-init modules-gs
  module load bcl-convert/4.2.7

  #
  # Infer the p5 barcode orientation from the flowcell RunParameters.xml file.
  #
  RCMP=`read_run_info.py ${illumina_run} RNA-seq`
echo \$RCMP >  /net/bbi/vol1/data/bge/bbi/tests/bbi-scirna-tests/rna3-065-a.rcmp/rcmp.txt
  #
  # Run bcl-convert.
  #
  $workflow.projectDir/bin/make_bclconvert_samplesheet.py -i $samplesheet_json -o samplesheet_bclconvert.csv -7 $p7_barcode_file -5 $p5_barcode_file --p5_rcmp \$RCMP
  bcl-convert \
    --bcl-input-directory $illumina_run \
    --output-directory fastqs_bclconvert \
    --sample-sheet samplesheet_bclconvert.csv
  """
}


