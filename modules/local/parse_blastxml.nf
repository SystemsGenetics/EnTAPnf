process PARSE_BLASTXML {
    label 'process_single'

    input:
    path blast_xml

    output:
    path "*.txt", emit: blast_txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    parse_blastxml.py --xml_file ${blast_xml} --out_file ${blast_xml}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_blastxml.py: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
