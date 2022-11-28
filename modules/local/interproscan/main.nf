process INTERPROSCAN {
    tag "$meta.id"
    label 'process_medium'

    // InterProScn has no conda package, a Galaxy singularity package, or a Quay,io docker image.
    // conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    // } else {
    //     container "quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0"
    // }
    container "annotater/interproscan:5.36-0.9"

    input:
    tuple val(meta), path(fasta)
    val seq_type

    output:
    tuple val(meta), path("*.{tsv,xml,gff,json,html,svg}"), emit: outfiles
    path "versions.yml"           , emit: version

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def stype = seq_type == 'pep' ? 'p' : 'n'
    """
    /usr/local/interproscan/interproscan.sh \\
        --input $fasta \\
        --cpu $task.cpus \\
        --output-file-base ${prefix} \\
        --seqtype $stype \\
        $args

    cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            \$(/usr/local/interproscan/interproscan.sh -version 2>&1 | head -n 1 | sed 's/^InterProScan version //')
    END_VERSIONS
    """
}
