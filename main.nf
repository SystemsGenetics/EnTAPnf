#!/usr/bin/env nextflow

/**
 * ========
 * GEMmaker
 * ========
 *
 * Authors:
 *  + John Hadish
 *  + Tyler Biggs
 *  + Stephen Ficklin
 *  + Ben Shealy
 *  + Connor Wytko
 *
 * Summary:
 *   A workflow for processing a large amount of RNA-seq data
 */



println """\
===================================
 FUNC-A   P I P E L I N E
===================================
General Information:
--------------------
  Profile(s):         ${workflow.profile}
  Container Engine:   ${workflow.containerEngine}

Input Parameters:
-----------------
  transcript (mRNA) file:     ${params.input.transcript_fasta}

Output Parameters:
------------------
  Output directory:           ${params.output.dir}

Execution Parameters:
---------------------
  Queue size:                 ${params.execution.queue_size}
  Number of threads:          ${params.execution.threads}
  Maximum retries:            ${params.execution.max_retries}
  Error strategy:             ${params.execution.error_strategy}

Software Parameters:
--------------------
  IPR core data:              ${params.software.interproscan.core_data_path}
  IPR panther data:           ${params.software.interproscan.panther_data_path}
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

  output:

  script:
    """
    /usr/local/interproscan/interproscan.sh -f TSV,XML,JSON,GFF3,HTML --goterms -i ${seq} --iprlookup --pathways --seqtype n 

    """
}
