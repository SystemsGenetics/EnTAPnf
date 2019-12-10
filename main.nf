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
  Profile(s):                 ${workflow.profile}
  Container Engine:           ${workflow.containerEngine}

Input Files:
-----------------
  Transcript FASTA file:      ${params.input.fasta_file}
  Transcript file type:       ${params.input.type}

Data Files:
-----------------
  InterProScan data:          ${params.data.interproscan}
  NCBI nr data:               ${params.data.nr}
  Uniprot SwissProt data:     ${params.data.sprot}
  OrthoDB data:               ${params.data.orthodb}
  STRING data:                ${params.data.string}

Output Parameters:
------------------
  Output directory:           ${params.output.dir}

"""
SEQS_FOR_IPRSCAN  = Channel.create()
SEQS_FOR_BLAST_NR = Channel.create()
SEQS_FOR_BLAST_SPROT = Channel.create()
SEQS_FOR_ORTHODB = Channel.create()
SEQS_FOR_BLAST_STRING = Channel.create()
SEQS_FOR_BLAST_TREMBL = Channel.create()

/**
 * Read in the transcript sequences. We will process them in small chunks
 * and injects those new sequence files as a tuple containing the file name
 * and the complete file into a set of channels for each major proecess.
 */
Channel.fromPath(params.input.fasta_file)
       .splitFasta(by: 10, file: true)
       .map {
         matches = it =~ /.*\/(.*)/
         [matches[0][1], it]
       }
       .separate(SEQS_FOR_IPRSCAN, SEQS_FOR_BLAST_NR, SEQS_FOR_BLAST_SPROT,SEQS_FOR_ORTHODB, SEQS_FOR_BLAST_STRING, SEQS_FOR_BLAST_TREMBL) { a -> [a, a, a, a, a, a]}


// Get the input sequence filename and put it in a value Channel so we
// can re-use it multiple times.
matches = params.input.fasta_file =~ /.*\/(.*)/
SEQUENCE_FILENAME = Channel.value(matches[0][1])

// Create a channel consisting of the array of OrthoDB level IDs
// to include in the analysis.
ORTHODB_LEVELS_LIST = Channel.from(params.steps.orthodb.levels)

// Create a channel indicating the type of blast to perform.
if (params.input.type == "nuc") {
  BLAST_TYPE = Channel.value('blastx')
  INTERPRO_TYPE = Channel.value('n')
}
else if (params.input.type == "pep") {
  BLAST_TYPE = Channel.value('blastp')
  INTERPRO_TYPE = Channel.value('p')
}
else {
  error "Error: the params.input.type setting must be either \"nuc\" for nucleotide or \"pep\" for peptide (i.e. protein sequence)."
}

// Make sure that if the SwissProt BLAST settings are good.
if (params.steps.dblast_sprot.enable == true) {
  // Make sure the data directory is present.
  data_file = file("${params.data.sprot}/uniprot_sprot.dmnd")
  if (data_file.isEmpty()) {
    error "Error: the Uniprot SwissProt Diamond index file cannot be found at ${params.data.sprot}/uniprot_sprot.dmnd. Please check the params.data.sprot setting and make sure the file is present in the specified directory."
  }
}

// Make sure that if the NR BLAST settings are good.
if (params.steps.dblast_nr.enable == true) {
  // Make sure the data directory is present.
  data_file = file("${params.data.nr}/nr.dmnd")
  if (data_file.isEmpty()) {
    error "Error: the NCBI nr Diamond index file cannot be found at ${params.data.nr}/nr.dmnd. Please check the params.data.nr setting and make sure the file is present in the specified directory."
  }
}

// Make sure that if the Trembl BLAST settings are good.
if (params.steps.dblast_trembl.enable == true) {
  // Make sure the data directory is present.
  data_file = file("${params.data.trembl}/uniprot_trembl.dmnd")
  if (data_file.isEmpty()) {
    error "Error: the Uniprot trembl Diamond index file cannot be found at ${params.data.trembl}/uniprot_trembl.dmnd. Please check the params.data.trembl setting and make sure the file is present in the specified directory."
  }
}

// Make sure that the String BLAST settings are good.
if (params.steps.string.enable == true) {
  // Make sure the data directory is present.
  data_file = file("${params.data.string}/protein.sequences.v11.0.dmnd")
  if (data_file.isEmpty()) {
    error "Error: the String Diamond index file cannot be found at ${params.data.string}/protein.sequences.v11.0.dmnd. Please check the params.data.string setting and make sure the file is present in the specified directory."
  }
}

// Make sure that if the OrthoDB database is specified the settings are good.
ORTHDB_TYPE  = Channel.create()
if (params.steps.orthodb.enable == true) {
  // Make sure the database name is valid.
  if (params.steps.orthodb.db != "plants" &&  params.steps.orthodb.db != "arthropoda" &&
      params.steps.orthodb.db != "verebrata" && params.steps.orthodb.db != "protozoa" &&
      params.steps.orthodb.db != "bacteria" && params.steps.orthodb.db != "fungi" &&
      params.steps.orthodb.db != "virdae") {
    error "Error: the params.steps.orthodb.db setting should be one of the following: \"plants\", \"arthropoda\", \"verebrata\", \"protozoa\", \"bacteria\", \"fungi\", or \"virdae\"."
  }
  ORTHDB_TYPE = Channel.value(params.steps.orthodb.db)

  // Make sure the data directory is present.
  data_dir = file("${params.data.orthodb}/${params.steps.orthodb.db}")
  if (data_dir.isEmpty()) {
    error "Error: the OrthoDB data directory cannot be found: ${params.data.orthodb}/${params.steps.orthodb.db}. Please check the params.data.orthodb setting."
  }
}

// Make sure that if the InterProScan settings are good.
if (params.steps.interproscan.enable == true) {
  // TODO: make sure the applicatin list is valid.

  // Make sure the data directory is present.
  data_dir = file("${params.data.interproscan}")
  if (data_dir.isEmpty()) {
    error "Error: the InterProScan data directory cannot be found: ${params.data.interproscan}. Please check the params.data.interproscan setting."
  }
}

/**
 * Retreives the list of species that should be searched for orthologs given
 * the levels provided by the user in the config file.
 */
process orthdb_level2species {
   label "orthdb_level2species"

   input:
     val level_id from ORTHODB_LEVELS_LIST

   output:
     stdout LEVELS_LIST_CSV

   when:
     params.steps.orthodb.enable == true

  script:
     """
     orthodb_level2species.py ${level_id} ${params.data.orthodb}/odb10v0_level2species.tab
     """
}
LEVELS_LIST_CSV.splitCsv().flatten().set{ORTHODB_SPECIES_LIST}


/**
 * Runs InterProScan on each sequence
 */
process interproscan {
  publishDir "${params.output.dir}/interpro/xml", mode: "link", pattern: "*.xml"
  label "interproscan"

  input:
    set val(seqname), file(seqfile) from SEQS_FOR_IPRSCAN
    val seq_type from INTERPRO_TYPE

  output:
    file "*.xml" into INTERPRO_XML
    file "*.tsv" into INTERPRO_TSV

  when:
    params.steps.interproscan.enable == true

  script:
    """
    # If this is kubernetes then hack a soft link for InterProScan's data
    if [ "${workflow.profile}" == "k8s" ]; then
      #EMPTY=""
      #if [ -n \${INTERPROSCAN_DATA_DIR+EMPTY} ]
      #then
      #    rm -fr /usr/local/interproscan/data
      #    ln -s \${INTERPROSCAN_DATA_DIR} /usr/local/interproscan/data
      #fi
      # Call InterProScan on a single sequence.
      echo "Hello!"
    fi
    /usr/local/interproscan/interproscan.sh \
      -f TSV,XML \
      --goterms \
      --input ${seqfile} \
      --iprlookup \
      --pathways \
      --seqtype ${seq_type} \
      --cpu ${task.cpus} \
      --output-dir . \
      --mode standalone \
      --applications ${params.steps.interproscan.applications}
    # Remove the temp directory created by InterProScan
    rm -rf ./temp
    """
}

/*
 * Wait until all Interpro jobs have finished and combine
 * all result files into a single lsit
 */
INTERPRO_TSV.collect().set{ INTERPRO_TSV_FILES }

/**
 * Combine InterProScan results.
 */
process interproscan_combine {
  publishDir "${params.output.dir}/interpro", mode: "link"
  label "interproscan_combine"

  input:
    file tsv_files from INTERPRO_TSV_FILES
    val sequence_filename from SEQUENCE_FILENAME

  output:
    file "${sequence_filename}.IPR_mappings.txt" into IPR_MAPPINGS
    file "${sequence_filename}.GO_mappings.txt" into GO_MAPPINGS
    file "${sequence_filename}.tsv" into INTEPRO_TSV_COMBINED

  script:
  """
    interpro_combine.py ${sequence_filename}
  """
}

/**
 * Runs Diamond blast against the NCBI non-redundant database.
 */
process dblast_nr {
  publishDir "${params.output.dir}/nr", mode: "link"
  label "diamond"
  memory "8 GB"

  input:
    val blast_type from BLAST_TYPE
    set val(seqname), file(seqfile) from SEQS_FOR_BLAST_NR

  output:
    file "*.xml" into BLAST_NR_XML

  when:
    params.steps.dblast_nr.enable == true

  script:
    """
    diamond ${blast_type} \
      --threads 1 \
      --query ${seqfile} \
      --db ${params.data.nr}/nr.dmnd \
      --out ${seqname}_vs_nr.${blast_type}.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}

/**
 * Runs Diamond blast against the Uniprot Trembl database.
 */
process dblast_trembl {
  publishDir "${params.output.dir}/trembl", mode: "link"
  label "diamond"
  memory "4 GB"

  input:
    val blast_type from BLAST_TYPE
    set val(seqname), file(seqfile) from SEQS_FOR_BLAST_TREMBL

  output:
    file "*.xml" into BLAST_TREMBL_XML

  when:
    params.steps.dblast_trembl.enable == true

  script:
    """
    diamond ${blast_type} \
      --threads 1 \
      --query ${seqfile} \
      --db ${params.data.trembl}/uniprot_trembl.dmnd \
      --out ${seqname}_vs_uniprot_trembl.${blast_type}.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}

/**
 * Runs Diamond blast against the SwissProt database.
 */
process dblast_sprot {
  publishDir "${params.output.dir}/sprot", mode: "link"
  label "diamond"

  input:
    val blast_type from BLAST_TYPE
    set val(seqname), file(seqfile) from SEQS_FOR_BLAST_SPROT

  output:
    file "*.xml"  into BLAST_SPROT_XML

  when:
    params.steps.dblast_sprot.enable == true

  script:
    """
    diamond ${blast_type} \
      --threads 1 \
      --query ${seqfile} \
      --db ${params.data.sprot}/uniprot_sprot.dmnd \
      --out ${seqname}_vs_uniprot_sprot.${blast_type}.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}

/**
 * Runs Diamond blast against the SwissProt database.
 */
process dblast_string {
  label "diamond"

  input:
    val blast_type from BLAST_TYPE
    set val(seqname), file(seqfile) from SEQS_FOR_BLAST_STRING

  output:
    file "*.xml"  into BLAST_STRING_XML

  when:
    params.steps.string.enable == true

  script:
    """
    diamond ${blast_type} \
      --threads 1 \
      --query ${seqfile} \
      --db ${params.data.string}/protein.sequences.v11.0.dmnd \
      --out ${seqname}_vs_string.${blast_type}.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}

/**
 * Runs Diamond blast against the OrthoDB species files.
 */
process dblast_orthodb {
  label "diamond"

  input:
    val blast_type from BLAST_TYPE
    set val(seqname), file(seqfile) from SEQS_FOR_ORTHODB

  output:
    file "*.xml" into ORTHODB_BLASTX_XML

  when:
    params.steps.orthodb.enable == true

  script:
    """
    diamond ${blast_type} \
      --threads 1 \
      --query ${seqfile} \
      --db ${params.data.orthodb}/odb10_all_og.dmnd \
      --out ${seqname}_vs_odb10_all_og.${blast_type}.xml \
      --evalue 1e-6 \
      --outfmt 5
    """
}

/**
 * Parses blast output against the orthdb blast databases and finds orthologs.
 */
process parse_dblast_orthodb {
   label "python3"

   input:
     file blast_xml from ORTHODB_BLASTX_XML

   output:
     file "*.txt" into ORTHODB_BLASTX_TXT

   when:
     params.steps.orthodb.enable == true

   script:
     """
     parse_blastxml.py --xml_file ${blast_xml} --out_file ${blast_xml}.txt
     """
}
