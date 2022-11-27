process PARSE_PROTEINPROTEIN {
    label 'process_single'

    container "annotater/python:3.7-0.9"

    input:
    path blast_xml
    val sequence_filename

    output:
    path "*.graphml", emit: graphml
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_proteinprotein.py \
        --xml ${blast_xml} \
        --db ${params.data.string}/protein \
        --species ${params.input.taxonomy_ID} \
        --out ${sequence_filename}.graphml

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_enzyme: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
