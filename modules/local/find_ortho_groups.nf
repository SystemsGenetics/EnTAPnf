process FIND_ORTHO_GROUPS {
    label 'process_single'

    container "annotater/python:3.7-0.9"

    input:
    path blast_xml
    val sequence_filename

    output:
    path "*.txt", emit: orths
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    identify_orthologous_groups.py \
        ${blast_xml} \
        ${params.data_orthodb} \
        ${sequence_filename}.orthodb_orthologs.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        identify_orthologous_groups.py: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
