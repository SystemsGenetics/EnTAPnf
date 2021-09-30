

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

    script:
    """
    interpro_combine.py ${sequence_filename}
    """
}
