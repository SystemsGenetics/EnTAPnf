process ENTAP_RUN {
    tag "$meta.id"
    label 'process_medium'

    // InterProScn has no conda package, a Galaxy singularity package, or a Quay,io docker image.
    // conda (params.enable_conda ? "bioconda::diamond=2.0.9" : null)
    // if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
    //     container 'https://depot.galaxyproject.org/singularity/diamond:2.0.9--hdcc8f71_0'
    // } else {
    //     container "quay.io/biocontainers/diamond:2.0.9--hdcc8f71_0"
    // }
    container "systemsgenetics/entap:flask"

    input:
    tuple val(meta), path(fasta)
    path (dbs)
    path (entap_config)
    val (seq_type)
    path (entap_outdir)

    output:
    tuple val(meta), path('*.tsv'), optional: true, emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_type = seq_type == 'pep' ? '--runP' : '--runN'
    def db_list = dbs.join(" -d ")
    """
    EnTAP $run_type \\
        -t $task.cpus  \\
        --ini $entap_config \\
        --input $fasta \\
        -d $db_list \\
        --out-dir $entap_outdir


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        diamond: \$(EnTAP --version 2>&1 | tail -n 1 | sed 's/^EnTAP  version: //')
    END_VERSIONS
    """
}
