process run_bclconvert {
  cache 'lenient'

  input:
    val samplesheet_json
    val illumina_run
    val p7_barcode_file
    val p5_barcode_file
    val p5_revcmp

  output:
    path "fastqs_bclconvert/*.fastq.gz"

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  #
  # Run bcl-convert.
  #
  $workflow.projectDir/bin/make_bclconvert_samplesheet.py -i $samplesheet_json -o samplesheet_bclconvert.csv -7 $p7_barcode_file -5 $p5_barcode_file --p5_rcmp $p5_revcmp

  bcl-convert \
    --bcl-input-directory $illumina_run \
    --output-directory fastqs_bclconvert \
    --sample-sheet samplesheet_bclconvert.csv
  """
}


