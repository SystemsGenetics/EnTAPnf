process INTERPROSCAN_COMBINE {
    label 'process_single'

    container "annotater/python:3.7-0.9"

    input:
    path tsv_files
    val sequence_filename

    output:
    path "${sequence_filename}.IPR_mappings.txt", emit: ipr_mappings
    path "${sequence_filename}.GO_mappings.txt", emit: go_mappings
    path "${sequence_filename}.tsv", emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    interpro_combine.py ${sequence_filename}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_enzyme: EnTAPnf ${workflow.workflow.manifest.version}
    END_VERSIONS
    """
}
