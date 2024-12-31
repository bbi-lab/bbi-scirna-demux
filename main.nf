import java.nio.file.Paths;

/*
** Run parameters.
**
** Moved to experiment.config.
**
params.samplesheet_json = '/net/bbi/vol1/data/bge/bbi/tests/bclconvert/RNA3-065-a/samplesheet.json'
params.illumina_run_dir = '/net/bbi/vol1/data/regression_tests/sciRNAseq/data_sources/illumina_runs/231121_NS500488_1448_AHJFNLBGXT'
params.output_dir = '/net/gs/vol1/home/bge/git/bbi-scirna-demux'
*/


/*
** 
*/
params.bin_dir = workflow.projectDir + '/bin'


/*
** Default barcode file paths.
*/
params.p7_barcode_file_default = "$workflow.projectDir/data/p7.txt"
params.p5_barcode_file_default = "$workflow.projectDir/data/p5.txt"
params.rt_barcode_file_default = "$workflow.projectDir/data/rt.txt"
params.ligation_barcode_file_default = "$workflow.projectDir/data/ligation.txt"


/*
** Set up channels.
*/
samplesheet_file = channel.value(params.samplesheet_json)
illumina_run_dir = channel.value(params.illumina_run_dir)
p7_barcode_file_default = channel.value(params.p7_barcode_file_default)
p5_barcode_file_default = channel.value(params.p5_barcode_file_default)
rt_barcode_file_default = channel.value(params.rt_barcode_file_default)
ligation_barcode_file_default = channel.value(params.ligation_barcode_file_default)


/*
** Import modules after defining params.* so that
** the parameters are accessible in the modules.
*/
include { run_bclconvert } from './modules/run_bclconvert.nf'
include { run_rna_rtlig_demux } from './modules/run_rna_rtlig_demux.nf'


/*
** Run pipeline.
*/
workflow {
  run_bclconvert(samplesheet_file, illumina_run_dir, p7_barcode_file_default, p5_barcode_file_default)
  run_bclconvert.out.flatMap{ make_pairwise_fastq_bclconvert(it) }.set{fastq_pairs}
  run_rna_rtlig_demux(fastq_pairs, samplesheet_file, rt_barcode_file_default, ligation_barcode_file_default)
}


/*
** Make tuples of fastq file pairs. The objects in fastqs_bclconvert
** are Unix paths. In order to substitute R2 for R1, convert to a
** string, replace, and convert to unix path.
*/
def make_pairwise_fastq_bclconvert(fastqs_bclconvert) {
  def fastq_pairs = []
  for(fastq_path in fastqs_bclconvert) {
    fastq_string = fastq_path.toString()
    if(fastq_string =~ '.+_R1_[0-9]{3}.fastq.gz') {
      if(fastq_string =~ '/Undetermined_.+') {
        continue
      }
      fastq_r1 = fastq_path
      fastq_r2 = Paths.get(fastq_string.replaceAll("_R1_", "_R2_"))
      tuple = new Tuple(fastq_r1, fastq_r2) 
      fastq_pairs.add(tuple)
    }
  }
  return(fastq_pairs)
}


