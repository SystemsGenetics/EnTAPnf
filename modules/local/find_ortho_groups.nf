
process FIND_ORTHO_GROUPS {
    publishDir "${params.outdir}/ortho_groups_orthodb",
        mode: params.publish_dir_mode

    container "annotater/python:3.7-0.9"

    input:
    file blast_xml
    val sequence_filename

    output:
    file "*.txt"

    script:
    """
    identify_orthologous_groups.py \
      ${blast_xml} \
      ${params.data_orthodb} \
      ${sequence_filename}.orthodb_orthologs.txt
    """
}
