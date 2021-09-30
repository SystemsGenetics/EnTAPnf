/*
========================================================================================
    VALIDATE INPUTS
========================================================================================
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowEntap.initialise(params, log)

// TODO nf-core: Add all file path parameters for the pipeline to the list below
// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input FASTA file not specified!' }

/*
========================================================================================
    CONFIG FILES
========================================================================================
*/


/*
========================================================================================
    IMPORT LOCAL MODULES/SUBWORKFLOWS
========================================================================================
*/

// Don't overwrite global params.modules, create a copy instead and use that within the main script.
def modules = params.modules.clone()

//
// MODULE: Local to the pipeline
//
include { GET_SOFTWARE_VERSIONS } from '../modules/local/get_software_versions' addParams( options: [publish_files : ['tsv':'']] )

/*
========================================================================================
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
========================================================================================
*/

//
// MODULE: Installed directly from nf-core/modules
//
def blastp_sprot_options = modules['diamond_blastp'].clone()
def blastx_sprot_options = modules['diamond_blastp'].clone()
blastp_sprot_options.publish_dir = 'blastp_vs_sprot'
blastx_sprot_options.publish_dir = 'blastx_vs_sprot'
include { DIAMOND_BLASTP as blastp_sprot } from '../modules/nf-core/modules/diamond/blastp/main.nf' addParams( options: blastp_sprot_options )
include { DIAMOND_BLASTX as blastx_sprot } from '../modules/nf-core/modules/diamond/blastx/main.nf' addParams( options: blastx_sprot_options )

def blastp_trembl_options = modules['diamond_blastp'].clone()
def blastx_trembl_options = modules['diamond_blastp'].clone()
blastp_trembl_options.publish_dir = 'blastp_vs_trembl'
blastx_trembl_options.publish_dir = 'blastx_vs_trembl'
include { DIAMOND_BLASTP as blastp_trembl } from '../modules/nf-core/modules/diamond/blastp/main.nf' addParams( options: blastp_trembl_options )
include { DIAMOND_BLASTX as blastx_trembl } from '../modules/nf-core/modules/diamond/blastx/main.nf' addParams( options: blastx_trembl_options )

def blastp_refseq_options = modules['diamond_blastp'].clone()
def blastx_refseq_options = modules['diamond_blastp'].clone()
blastp_refseq_options.publish_dir = 'blastp_vs_refseq'
blastx_refseq_options.publish_dir = 'blastx_vs_refseq'
include { DIAMOND_BLASTP as blastp_refseq } from '../modules/nf-core/modules/diamond/blastp/main.nf' addParams( options: blastp_refseq_options )
include { DIAMOND_BLASTX as blastx_refseq } from '../modules/nf-core/modules/diamond/blastx/main.nf' addParams( options: blastx_refseq_options )

def blastp_nr_options = modules['diamond_blastp'].clone()
def blastx_nr_options = modules['diamond_blastp'].clone()
blastp_nr_options.publish_dir = 'blastp_vs_nr'
blastx_nr_options.publish_dir = 'blastx_vs_nr'
include { DIAMOND_BLASTP as blastp_nr } from '../modules/nf-core/modules/diamond/blastp/main.nf' addParams( options: blastp_nr_options )
include { DIAMOND_BLASTX as blastx_nr } from '../modules/nf-core/modules/diamond/blastx/main.nf' addParams( options: blastx_nr_options )

def blastp_orthodb_options = modules['diamond_blastp'].clone()
def blastx_orthodb_options = modules['diamond_blastp'].clone()
blastp_orthodb_options.publish_dir = 'blastp_vs_orthodb'
blastx_orthodb_options.publish_dir = 'blastx_vs_orthodb'
include { DIAMOND_BLASTP as blastp_orthodb } from '../modules/nf-core/modules/diamond/blastp/main.nf' addParams( options: blastp_orthodb_options )
include { DIAMOND_BLASTX as blastx_orthodb } from '../modules/nf-core/modules/diamond/blastx/main.nf' addParams( options: blastx_orthodb_options )

def blastp_string_options = modules['diamond_blastp'].clone()
def blastx_string_options = modules['diamond_blastp'].clone()
blastp_string_options.publish_dir = 'blastp_vs_string'
blastx_string_options.publish_dir = 'blastx_vs_string'
include { DIAMOND_BLASTP as blastp_string } from '../modules/nf-core/modules/diamond/blastp/main.nf' addParams( options: blastp_string_options )
include { DIAMOND_BLASTX as blastx_string } from '../modules/nf-core/modules/diamond/blastx/main.nf' addParams( options: blastx_string_options )

def interproscan_nuc_options = modules['interproscan'].clone()
def interproscan_pep_options = modules['interproscan'].clone()
interproscan_nuc_options.args = interproscan_nuc_options.args + "--seqtype n "
interproscan_nuc_options.args = interproscan_nuc_options.args.replace("/--applications '.*?'", "---applications ${params.ipr_apps}")
interproscan_pep_options.args = interproscan_pep_options.args.replace("/--applications '.*?'", "---applications ${params.ipr_apps}")
include { INTERPROSCAN as interproscan_nuc } from '../modules/local/modules/interproscan/main.nf' addParams( options: interproscan_nuc_options )
include { INTERPROSCAN as interproscan_pep } from '../modules/local/modules/interproscan/main.nf' addParams( options: interproscan_pep_options )

include { INTERPROSCAN_COMBINE as interproscan_combine } from '../modules/local/interproscan_combine.nf'
include { FIND_EC_NUMBERS as find_ec_numbers } from '../modules/local/find_EC_numbers.nf'
include { FIND_ORTHO_GROUPS as find_ortho_groups } from '../modules/local/find_ortho_groups.nf'

/*
========================================================================================
    RUN MAIN WORKFLOW
========================================================================================
*/


workflow ENTAP {

    ch_software_versions = Channel.empty()


    //
    // MODULE: Pipeline reporting
    //
    ch_software_versions
        .map { it -> if (it) [ it.baseName, it ] }
        .groupTuple()
        .map { it[1][0] }
        .flatten()
        .collect()
        .set { ch_software_versions }

    //
    // Split the FASTA file into groups of size 10.
    //
    sequence_filename =  params.input.replaceFirst(/^.*\/(.*)$/, '$1')
    ch_split_seqs = Channel
        .fromPath(params.input)
        .splitFasta(by: params.batch_size, file: true)
        // Add an ID for each split file which is essential the file name
        .map {
            [[id: it.getFileName()], it]
        }
        // Because we're using Diamond against many different databases
        // we need to add a quantifier to the ID to indicate the database
        // that will be used. Otherwise the output name won't indicate the
        // database and we'll have conflicts.  The following creates a
        // sub channels for each tool.
        .multiMap { it ->
            sprot: [[id: it[0].id  + '_vs_sprot'], it[1]]
            nr: [[id: it[0].id  + '_vs_nr'], it[1]]
            refseq: [[id: it[0].id  + '_vs_refseq'], it[1]]
            orthodb: [[id: it[0].id  + '_vs_orthodb'], it[1]]
            string: [[id: it[0].id  + '_vs_string'], it[1]]
            trembl: [[id: it[0].id  + '_vs_trembl'], it[1]]
            ipr: it
        }
    //ch_split_seqs.sprot.view()

    GET_SOFTWARE_VERSIONS (
        ch_software_versions.map { it }.collect()
    )

    //
    // BLAST sequences against ExPASy SwissProt using Diamond.
    //
    if (params.data_sprot) {
        db = [ file(params.data_sprot, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_sprot(ch_split_seqs.sprot, db)
            blastp_sprot.out.txt
                .map { it[1] }
                .set { blastp_sprot_txt }
            find_ec_numbers(blastp_sprot_txt, sequence_filename)
        }
        if (params.seq_type == 'nuc') {
            blastx_sprot(ch_split_seqs.sprot, db)
            blastx_sprot.out.txt
                .map { it[1] }
                .set { blastx_sprot_txt }
            find_ec_numbers(blastx_sprot_txt, sequence_filename)
        }
    }


    //
    // BLAST sequences against ExPASy tembl using Diamond.
    //
    if (params.data_trembl) {
        db = [ file(params.data_trembl, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_data_trembl(ch_split_seqs.data_trembl, db)
        }
        if (params.seq_type == 'nuc') {
            blastx_data_trembl(ch_split_seqs.data_trembl, db)
        }
    }

    //
    // BLAST sequences against NCBI nr using Diamond.
    //
    if (params.data_nr) {
        db = [ file(params.data_sprot, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_nr(ch_split_seqs.sprot, db)
        }
        if (params.seq_type == 'nuc') {
            blastx_nr(ch_split_seqs.sprot, db)
        }
    }

    //
    // BLAST sequences against NCBI refseq using Diamond.
    //
    if (params.data_refseq) {
        db = [ file(params.data_sprot, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_refseq(ch_split_seqs.refseq, db)
        }
        if (params.seq_type == 'nuc') {
            blastx_refseq(ch_split_seqs.refseq, db)
        }
    }

    //
    // BLAST sequences against OrthoDB using Diamond.
    //
    if (params.data_orthodb) {
        db = [ file(params.data_orthodb, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_orthodb(ch_split_seqs.orthodb, db)
            blastp_orthodb.out.txt
                .map { it[1] }
                .set { blastp_orthodb_txt }
            find_ortho_groups(blastp_orthodb_txt, sequence_filename)
        }
        if (params.seq_type == 'nuc') {
            blastx_orthodb(ch_split_seqs.orthodb, db)
            blastx_orthodb.out.txt
                .map { it[1] }
                .set { blastx_orthodb_txt }
            find_ortho_groups(blastx_orthodb_txt, sequence_filename)
        }
    }

    //
    // BLAST sequences against String using Diamond.
    //
    if (params.data_string) {
        db = [ file(params.data_string, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_string(ch_split_seqs.string, db)
        }
        if (params.seq_type == 'nuc') {
            blastx_string(ch_split_seqs.string, db)
        }
    }
    //
    // InterProScan
    //
    if (params.data_ipr) {

        if (params.seq_type == 'pep') {
            interproscan_pep(ch_split_seqs.ipr)
            interproscan_pep.out.outfiles
                .map { it[1] }
                .flatten()
                .branch {
                    tsv: it.getFileName().toString().endsWith(".tsv")
                    xml: it.getFileName().toString().endsWith(".xml")
                }
                .set { interproscan_pep_out }
            interproscan_combine(interproscan_pep_out.tsv, sequence_filename)
        }
        if (params.seq_type == 'nuc') {
            interproscan_nuc(ch_split_seqs.ipr, fromPath(params.input))
            interproscan_nuc.out.outfiles
                .map { it[1] }
                .flatten()
                .branch {
                    tsv: it.getFileName().toString().endsWith(".tsv")
                    xml: it.getFileName().toString().endsWith(".xml")
                }
                .set { interproscan_pep_out }
            interproscan_combine(interproscan_nuc_out.tsv, sequence_filename)
        }


    }
}

/*
========================================================================================
    COMPLETION EMAIL AND SUMMARY
========================================================================================
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
}

/*
========================================================================================
    THE END
========================================================================================
*/
