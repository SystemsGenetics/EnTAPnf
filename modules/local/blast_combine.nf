process BLAST_COMBINE {
    label 'process_single'

    container "annotater/python:3.7-0.9"

    input:
    path outfiles
    val blast_type
    val sequence_filename
    val db_name

    output:
    path ("${blast_type}_${sequence_filename}_${db_name}.out"), emit: out
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    combine_blast.sh "${blast_type}_${sequence_filename}_${db_name}.out"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        interpro_combine.py: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
