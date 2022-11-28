process FIND_EC_NUMBERS {
    tag "$meta.id"
    label 'process_single'

    container "systemsgenetics/entap:flask"

    input:
    tuple val(meta), path(blast_xml)
    path (enzyme_dat)

    output:
    tuple val(meta), path ("*.txt"), emit: ecout
    path ("versions.yml"), emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    parse_enzyme.py \
        --xml $blast_xml \
        --enzyme $enzyme_dat \
        --out ${prefix}.ECnumbers.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        parse_enzyme.py: EnTAPnf ${workflow.manifest.version}
    END_VERSIONS
    """
}
