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

/**
 * Read in the transcript sequences. We will process them one at a time.
 */
if (params.input.transcript_fasta == "none") {
  Channel.empty().set { TRANSCRIPT_SEQS_FOR_IPRSCAN }
}
else {
  Channel.fromPath(params.input.transcript_fasta)
    .splitFasta(by: 1, file: true)
    .set { TRANSCRIPT_SEQS_FOR_IPRSCAN }
}


/**
 * Runs InterProScan on each sequence
 */
process interproscan {
  publishDir params.output.dir, mode: "symlink"
  label "interproscan"

  input:
    file seq from TRANSCRIPT_SEQS_FOR_IPRSCAN

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
