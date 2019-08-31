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

Input Parameters:
-----------------
  Transcript (mRNA) file:     ${params.input.transcript_fasta}
  InterProScan data:          ${params.input.interproscan}
  Panther data:               ${params.input.panther}
  NCBI nr data:               ${params.input.nr}
  Uniprot SwissProt data:     ${params.input.uniprot_sprot}

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


/**
 * Runs InterProScan on each sequence
 */
process interproscan {
  publishDir params.output.dir, mode: "symlink"
  label "interproscan"

  input:
    file seq from SEQS_FOR_IPRSCAN

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
      --applications ${params.software.interproscan.applications}
    # Remove the temp directory created by InterProScan
    rm -rf ./temp
    """
}

/**
 * Runs NCBI blastx against the NCBI non-redundant database.
 */
process blastx_nr {
  publishDir params.output.dir, mode: "symlink"
  label "ncbi_blast"

  input:
    file seq from SEQS_FOR_BLASTX_NR

  script:
    """
    /usr/local/ncbi-blast/bin/blastx \
      -query ${seq} \
      -db /annotater/nr/nr \
      -out ${seq}.blastx.xml \
      -evalue 1e-6 \
      -outfmt 13
    """
}

/**
 * Runs NCBI blastx against the SwissProt database.
 */
process blastx_sprot {
  publishDir params.output.dir, mode: "symlink"
  label "ncbi_blast"

  input:
    file seq from SEQS_FOR_BLASTX_SPROT

  script:
    """
    /usr/local/ncbi-blast/bin/blastx \
      -query ${seq} \
      -db /annotater/uniprot_sprot/uniprot_sprot \
      -out ${seq}.blastx.xml \
      -evalue 1e-6 \
      -outfmt 13
    """
}
