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
    path (entap_db)
    path (eggnog_db)
    path (data_eggnog)
    path (blast_results)
    path (interproscan_results)

    output:
    tuple val(meta), path('*.tsv'), optional: true, emit: tsv
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def run_type = seq_type == 'pep' ? '--runP' : '--runN'
    def db_list = dbs.join(" -d \$PWD/")
    """
    # Update the config file so it can find the eggnong diamond file
    cp $entap_config new.$entap_config
    CWD=`echo \$PWD | perl -p -e 's/\\\\//\\\\\\\\\\\\//g'`
    echo "s/eggnog-dmnd=eggnog_proteins.dmnd/eggnog-dmnd=\$CWD\\\\/$data_eggnog/"
    perl -pi -e "s/eggnog-dmnd=eggnog_proteins.dmnd/eggnog-dmnd=\$CWD\\\\/$data_eggnog/" new.$entap_config
    perl -pi -e "s/eggnog-sql=eggnog.db/eggnog-sql=\$CWD\\\\/$eggnog_db/" new.$entap_config

    # Link the blast files to the directory EnTAP expects so it doesn't
    # rerun those.
    mkdir -p \$PWD/outfiles/similarity_search/DIAMOND
    files=`ls blast*.out`
    for f in \$files; do
        ln -s \$PWD/\$f \$PWD/outfiles/similarity_search/DIAMOND/\$f
    done;

    # Link the InterProScan file to the diretory EnTAP expects
    ln -s \$PWD/$interproscan_results \$PWD/ontology/InterProScan/interpro_results.tsv

    # Run EnTAP
    EnTAP $run_type \\
        -t $task.cpus \\
        --ini new.$entap_config \\
        --input $fasta \\
        -d \$PWD/$db_list \\
        --out-dir \$PWD/outfiles


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        EnTAP: \$(EnTAP --version 2>&1 | tail -n 1 | sed 's/^EnTAP  version: //')
    END_VERSIONS
    """
}
