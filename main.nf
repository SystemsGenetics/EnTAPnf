#!/usr/bin/env nextflow

/**
 * ========
 * AnnoTater
 * ========
 *
 * Authors:
 *  + Stephen Ficklin
 *
 * Summary:
 *   A workflow for annotating Eukaryotic transcript sequences from whole genome
 *   or de novo transcriptome assemblies.
 */



println """\

General Information:
--------------------
  Profile(s):         ${workflow.profile}
  Container Engine:   ${workflow.containerEngine}

Input Files:
-----------------
  Transcript (mRNA) file:     ${params.input.transcript_fasta}

Data Files:
-----------------
  InterProScan data:          ${params.data.interproscan}
  Panther data:               ${params.data.panther}
  NCBI nr data:               ${params.data.nr}
  Uniprot SwissProt data:     ${params.data.uniprot_sprot}

Output Parameters:
------------------
  Output directory:           ${params.output.dir}

"""
SEQS_FOR_IPRSCAN  = Channel.create()
SEQS_FOR_BLASTX_NR = Channel.create()
SEQS_FOR_BLASTX_SPROT = Channel.create()

/**
 * Read in the transcript sequences. We will process them one at a time.
 */
Channel.fromPath(params.input.transcript_fasta)
       .splitFasta(by: 1, file: true)
       .separate(SEQS_FOR_IPRSCAN, SEQS_FOR_BLASTX_NR, SEQS_FOR_BLASTX_SPROT) { a -> [a, a, a]}


process orthodb_index {
  label "diamond_makedb"
  cpus = 2

  output:
    file "*.dmnd" into ORTHODB_INDEXES

  when:
    params.steps.orthodb.enable == true

  script:
    if (params.data.orthodb.dbs.plants == true)
      """
      for species in ${params.data.orthodb.species.join(" ")}; do
        diamond makedb \
          --threads 2 \
          --in /annotater/orthodb/plants/Rawdata/\${species}_0.fs \
          --db \${species}
      done
      """
    else
      """
        echo "A database for OrthoDB has not been selected"
        echo ${params.data.orthodb.dbs.plants}
        exit 1
      """
}
/**
 * Prepares the Diamond indexes for the Uniprot Sprot database.
 */
process uniprot_sprot_index {
  label "diamond_makedb"
  cpus = 2

  output:
    file "*.dmnd" into SPROT_INDEX

  when:
    params.steps.dblastx_sprot.enable == true

  script:
  """
    diamond makedb \
      --threads 2 \
      --in /annotater/uniprot_sprot/uniprot_sprot.fasta \
      --db uniprot_sprot
  """
}
/**
 * Prepares the Diamond indexes for the NCBI nr database.
 */
process nr_index {
  label "diamond_makedb"
  cpus = 2
  memory = "6 GB"

  output:
    file "*.dmnd" into NR_INDEX

  when:
    params.steps.dblastx_nr.enable == true

  script:
  """
    diamond makedb \
      --threads 2 \
      --index-chunks 100 \
      --in /annotater/nr/nr \
      --db nr
  """
}

/**
 * Runs InterProScan on each sequence
 */
process interproscan {
  label "interproscan"

  input:
    file seq from SEQS_FOR_IPRSCAN

  output:
    file "*.json" into INTERPRO_JSON

  when:
    params.steps.interproscan.enable == true

  script:
    """
    # Call InterProScan on a single sequence.
    /usr/local/interproscan/interproscan.sh \
      -f JSON \
      --goterms \
      --input ${seq} \
      --iprlookup \
      --pathways \
      --seqtype n \
      --cpu ${task.cpus} \
      --output-dir . \
      --mode standalone \
      --applications ${params.steps.interproscan.applications}
    # Remove the temp directory created by InterProScan
    rm -rf ./temp
    """
}

/**
 * Runs blastx against the NCBI non-redundant database.
 */
process dblastx_nr {
  label "diamond"

  input:
    file seq from SEQS_FOR_BLASTX_NR
    file index from NR_INDEX

  output:
    file "*_vs_nr.dblastx.xml" into BLASTX_NR_XML

  when:
    params.steps.dblastx_nr.enable == true

  script:
    """
    diamond blastx \
      --threads 1 \
      --query ${seq} \
      --db nr \
      --out ${seq}_vs_nr.dblastx.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}

/**
 * Runs blastx against the SwissProt database.
 */
process dblastx_sprot {
  label "diamond"

  input:
    file seq from SEQS_FOR_BLASTX_SPROT
    file index from SPROT_INDEX

  output:
    file "*_vs_uniprot_sprot.blastx_1.xml"  into BLASTX_SPROT_XML

  when:
    params.steps.dblastx_sprot.enable == true

  script:
    """
    diamond blastx \
      --threads 1 \
      --query ${seq} \
      --db uniprot_sprot \
      --out ${seq}_vs_uniprot_sprot.blastx.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}
