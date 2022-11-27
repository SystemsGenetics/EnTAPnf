process FIND_EC_NUMBERS {
    label 'process_single'

    container "annotater/python:3.7-0.9"

    input:
    path blast_xml
    val sequence_filename

    output:
    path "*.txt", emit: ecout
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_enzyme.py \
        --xml ${blast_xml} \
        --enzyme ${params.data_sprot}/enzyme.dat \
        --out ${sequence_filename}.ECnumbers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_enzyme: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
