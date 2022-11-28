process ENTAP_CONFIG {
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
    path (entap_config)

    output:
    path('outfiles/databases/entap_database.db'), emit: entap_db
    path('outfiles/databases/eggnog.db'), emit: eggnog_db
    path('outfiles/bin/eggnog_proteins.dmnd'), emit: data_eggnog
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    EnTAP --config \\
        -t $task.cpus  \\
        --ini $entap_config \\
        --out-dir ./outfiles


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EnTAP: \$(EnTAP --version 2>&1 | tail -n 1 | sed 's/^EnTAP  version: //')
    END_VERSIONS
    """
}
