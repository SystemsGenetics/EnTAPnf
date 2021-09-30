// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)


process INTERPROSCAN_COMBINE {
    publishDir "${params.outdir}/interproscan_combined",
        mode: params.publish_dir_mode

    container "annotater/python:3.7-0.9"

    input:
    file tsv_files
    val sequence_filename

    output:
    file "${sequence_filename}.IPR_mappings.txt"
    file "${sequence_filename}.GO_mappings.txt"
    file "${sequence_filename}.tsv"
    path "versions.yml"           , emit: version

    script:
    """
    interpro_combine.py ${sequence_filename}

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: 'EnTAP-nf internal script')
    END_VERSIONS
    """
}
