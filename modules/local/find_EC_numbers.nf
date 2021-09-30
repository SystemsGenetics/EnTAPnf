
process FIND_EC_NUMBERS {
    publishDir "${params.outdir}/EC_numbers_sprot",
        mode: params.publish_dir_mode

    container "annotater/python:3.7-0.9"

    input:
    file blast_xml
    val sequence_filename

    output:
    file "*.txt"

    script:
    """
    parse_enzyme.py \
        --xml ${blast_xml} \
        --enzyme ${params.data_sprot}/enzyme.dat \
        --out ${sequence_filename}.ECnumbers.txt
    """
}
