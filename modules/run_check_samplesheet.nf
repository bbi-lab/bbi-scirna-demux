process run_check_samplesheet {
  cache 'lenient'

  input:
    val samplesheet_json

  script:
  """
  # bash watch for errors
  set -ueo pipefail

  $workflow.projectDir/bin/check_samplesheet_json.py -i ${samplesheet_json}
 
  """
}


