// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process INTERPROSCAN {
    tag "$meta.id"
    label 'process_medium'
    label 'interproscan'
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:getSoftwareName(task.process), meta:meta, publish_by_meta:['id']) }

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

    output:
    tuple val(meta), path("*.{tsv,xml,gff,json,html,svg}"), emit: outfiles
    path "versions.yml"           , emit: version

    script:
    def software = getSoftwareName(task.process)
    def prefix   = options.suffix ? "${meta.id}${options.suffix}" : "${meta.id}_interpro"
    """
    /usr/local/interproscan/interproscan.sh \\
          --input $fasta \\
          --cpu $task.cpus \\
          --output-file-base ${prefix} \\
          $options.args

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        ${getSoftwareName(task.process)}: \$(/usr/local/interproscan/interproscan.sh -version 2>&1 | head -n 1 | sed 's/^InterProScan version //')
    END_VERSIONS
    """
}
