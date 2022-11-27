/**
 * Parses blast output against protein protein interaction.
 */
process PARSE_PROTEINPROTEIN {
    publishDir "${params.output.dir}/proteinprotein/", mode: "link"

    container "annotater/python:3.7-0.9"

    input:
        file blast_xml from BLAST_STRING_XML
        val sequence_filename from SEQUENCE_FILENAME

    output:
        file "*.graphml" into PROTEINPROTEIN_BLASTX_TXT

    script:
    """
    parse_proteinprotein.py \
        --xml ${blast_xml} \
        --db ${params.data.string}/protein \
        --species ${params.input.taxonomy_ID} \
        --out ${sequence_filename}.graphml
    """
}
