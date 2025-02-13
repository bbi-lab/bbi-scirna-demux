process run_bclconvert {
  cache 'lenient'

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
  # I believe that bcl-convert manages the P5 orientation.
  #
  RCMP='False'

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


