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

  $workflow.projectDir/bin/make_bclconvert_samplesheet.py -i $samplesheet_json -o samplesheet_bclconvert.csv -7 $p7_barcode_file -5 $p5_barcode_file
  bcl-convert \
    --bcl-input-directory $illumina_run \
    --output-directory fastqs_bclconvert \
    --sample-sheet samplesheet_bclconvert.csv
  """
}


