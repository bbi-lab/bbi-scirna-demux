def demux_out = params.output_dir + '/demux_out'
def demux_log = params.output_dir + '/demux_log'


process run_rna_rtlig_demux {
  cache 'lenient'

  publishDir path: "${demux_out}", pattern: "*.bam", mode: 'copy'
  publishDir path: "${demux_log}", pattern: "*.json", mode: 'copy'

  input:
  tuple path(fastq_read1), path(fastq_read2)
  val(samplesheet_file)
  val(rt_file)
  val(ligation_file)

  output:
  path("*.bam")
  path("*.json")

  script:
  """
  $workflow.projectDir/bin/rna_rtlig_demux -1 $fastq_read1 \
                                           -2 $fastq_read2 \
                                           -s $samplesheet_file \
                                           -r $rt_file \
                                           -l $ligation_file \
                                           -f bam \
                                           --ncpu 3
  """
}


