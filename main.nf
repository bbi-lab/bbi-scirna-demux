import java.nio.file.Paths;

/*
** Run parameters.
*/
params.p5_revcmp = 'False'


/*
** 
*/
params.bin_dir = workflow.projectDir + '/bin'
params.recipe = 'std_1'      // barcode 'recipe'. 'std_1' is the default, which is for non-MegaSci. 'std_2' is for MegaSci.


/*
** Default barcode file paths.
*/
params.p7_barcode_file_default = "$workflow.projectDir/data/p7.txt"
params.p5_barcode_file_default = "$workflow.projectDir/data/p5.txt"
params.rt_barcode_file_default = "$workflow.projectDir/data/rt.txt"
params.ligation_barcode_file_default = "$workflow.projectDir/data/ligation.txt"


/*
** Import modules after defining params.* so that
** the parameters are accessible in the modules.
*/
include {run_check_samplesheet} from './modules/run_check_samplesheet.nf'
include { run_bclconvert } from './modules/run_bclconvert.nf'
include { run_rna_rtlig_demux } from './modules/run_rna_rtlig_demux.nf'


/*
** Run pipeline.
*/
workflow {
  def samplesheet_file = channel.value(params.samplesheet_json)
  def illumina_run_dir = channel.value(params.illumina_run_dir)
  def p7_barcode_file_default = channel.value(params.p7_barcode_file_default)
  def p5_barcode_file_default = channel.value(params.p5_barcode_file_default)
  def rt_barcode_file_default = channel.value(params.rt_barcode_file_default)
  def ligation_barcode_file_default = channel.value(params.ligation_barcode_file_default)

  run_check_samplesheet(samplesheet_file)
  run_bclconvert(samplesheet_file, illumina_run_dir, p7_barcode_file_default, p5_barcode_file_default, params.p5_revcmp)
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
  for(def fastq_path in fastqs_bclconvert) {
    def fastq_string = fastq_path.toString()
    if(fastq_string =~ '.+_R1_[0-9]{3}.fastq.gz') {
      if(fastq_string =~ '/Undetermined_.+') {
        continue
      }
      def fastq_r1 = fastq_path
      def fastq_r2 = Paths.get(fastq_string.replaceAll("_R1_", "_R2_"))
      def tuple = new Tuple(fastq_r1, fastq_r2) 
      fastq_pairs.add(tuple)
    }
  }
  return(fastq_pairs)
}


