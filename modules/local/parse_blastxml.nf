process PARSE_BLASTXML {
    tag "$meta.id"
    label 'process_single'

    container "systemsgenetics/entap:flask"

    input:
    tuple val(meta), path(blast_xml)

    output:
    tuple val(meta), path ("*.txt"), emit: blast_txt
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    outfile=`basename ${blast_xml} .xml`
    parse_blastxml.py --xml_file ${blast_xml} --out_file \${outfile}.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_blastxml.py: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
