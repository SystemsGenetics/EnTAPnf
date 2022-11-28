/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE INPUTS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def summary_params = NfcoreSchema.paramsSummaryMap(workflow, params)

// Validate input parameters
WorkflowEntapnf.initialise(params, log)

// Check input path parameters to see if they exist
def checkPathParamList = [ params.input ]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

// Check mandatory parameters
if (params.input) { ch_input = file(params.input) } else { exit 1, 'Input FASTA file not specified!' }

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { INTERPROSCAN as interproscan_nuc } from '../modules/local/interproscan/main.nf'
include { INTERPROSCAN as interproscan_pep } from '../modules/local/interproscan/main.nf'
include { INTERPROSCAN_COMBINE as interproscan_combine } from '../modules/local/interproscan_combine.nf'
include { FIND_EC_NUMBERS as find_ec_numbers } from '../modules/local/find_EC_numbers.nf'
include { FIND_ORTHO_GROUPS as find_ortho_groups } from '../modules/local/find_ortho_groups.nf'
include { ENTAP_CONFIG as entap_config } from '../modules/local/entap/entap_config/main.nf'
include { ENTAP_RUN as entap_run } from '../modules/local/entap/entap_run/main.nf'

include { PARSE_BLASTXML as parse_blastp_sprot} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastx_sprot} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastp_trembl} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastx_trembl} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastp_refseq} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastx_refseq} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastp_nr} from '../modules/local/parse_blastxml.nf'
include { PARSE_BLASTXML as parse_blastx_nr} from '../modules/local/parse_blastxml.nf'

include { BLAST_COMBINE as combine_blastp_sprot} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastx_sprot} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastp_trembl} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastx_trembl} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastp_refseq} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastx_refseq} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastp_nr} from '../modules/local/blast_combine.nf'
include { BLAST_COMBINE as combine_blastx_nr} from '../modules/local/blast_combine.nf'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'

//
// MODULE: Installed directly from nf-core/modules
//
include { DIAMOND_BLASTP as blastp_sprot } from '../modules/nf-core/diamond/blastp/main.nf'
include { DIAMOND_BLASTX as blastx_sprot } from '../modules/nf-core/diamond/blastx/main.nf'
include { DIAMOND_BLASTP as blastp_trembl } from '../modules/nf-core/diamond/blastp/main.nf'
include { DIAMOND_BLASTX as blastx_trembl } from '../modules/nf-core/diamond/blastx/main.nf'
include { DIAMOND_BLASTP as blastp_refseq } from '../modules/nf-core/diamond/blastp/main.nf'
include { DIAMOND_BLASTX as blastx_refseq } from '../modules/nf-core/diamond/blastx/main.nf'
include { DIAMOND_BLASTP as blastp_nr } from '../modules/nf-core/diamond/blastp/main.nf'
include { DIAMOND_BLASTX as blastx_nr } from '../modules/nf-core/diamond/blastx/main.nf'
include { DIAMOND_BLASTP as blastp_orthodb } from '../modules/nf-core/diamond/blastp/main.nf'
include { DIAMOND_BLASTX as blastx_orthodb } from '../modules/nf-core/diamond/blastx/main.nf'
include { DIAMOND_BLASTP as blastp_string } from '../modules/nf-core/diamond/blastp/main.nf'
include { DIAMOND_BLASTX as blastx_string } from '../modules/nf-core/diamond/blastx/main.nf'


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow ENTAPNF {

    ch_versions = Channel.empty()

    //
    // Split the FASTA file into groups of size 10.
    //
    sequence_filename = params.input.replaceFirst(/^.*\/(.*)\..+$/, '$1')
    entap_sequence_filename = params.input.replaceFirst(/^.*\/([^\.]+).*$/, '$1') + "_final"
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

    // Houses the list of diamond database files for EnTAP
    entap_dbs = []

    //
    // BLAST sequences against ExPASy SwissProt using Diamond.
    //
    if (params.data_sprot) {
        db = [ file(params.data_sprot + '/uniprot_sprot.dmnd', checkIfExists: true) ]
        entap_dbs.add(db[0])
        if (params.seq_type == 'pep') {
            blastp_sprot(ch_split_seqs.sprot, db, 'xml', '')
            if (params.enzyme_dat) {
                find_ec_numbers(blastp_sprot.out.xml, params.enzyme_dat)
            }
            parse_blastp_sprot(blastp_sprot.out.xml)
            parse_blastp_sprot.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastp_sprot_txt }
            combine_blastp_sprot(blastp_sprot_txt, 'blastp', entap_sequence_filename, 'uniprot_sprot')
        }
        if (params.seq_type == 'nuc') {
            blastx_sprot(ch_split_seqs.sprot, db, 'xml', '')
            if (params.enzyme_dat) {
                find_ec_numbers(blastx_sprot.out.xml, params.enzyme_dat)
            }
            parse_blastx_sprot(blastx_sprot.out.xml)
            parse_blastx_sprot.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastx_sprot_txt }
            combine_blastx_sprot(blastx_sprot_txt, 'blastx', entap_sequence_filename, 'uniprot_sprot')
        }
    }


    //
    // BLAST sequences against ExPASy tembl using Diamond.
    //
    if (params.data_trembl) {
        db = [ file(params.data_trembl + "/uniprot_trembl.dmnd", checkIfExists: true) ]
        entap_dbs.add(db[0])
        if (params.seq_type == 'pep') {
            blastp_data_trembl(ch_split_seqs.data_trembl, db, 'xml', '')
            parse_blastp_trembl(blastp_trembl.out.xml)
            parse_blastp_trembl.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastp_trembl_txt }
            combine_blastp_trembl(blastp_ntrembl_txt, 'blastp', entap_sequence_filename, 'uniprot_trembl')
        }
        if (params.seq_type == 'nuc') {
            blastx_data_trembl(ch_split_seqs.data_trembl, db, 'xml', '')
            parse_blastx_trembl(blastx_trembl.out.xml)
            parse_blastx_trembl.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastx_trembl_txt }
            combine_blastx_trembl(blastx_ntrembl_txt, 'blastx', entap_sequence_filename, 'uniprot_trembl')
        }
    }

    //
    // BLAST sequences against NCBI nr using Diamond.
    //
    if (params.data_nr) {
        db = [ file(params.data_nr + "/nr.dmnd", checkIfExists: true) ]
        entap_dbs.add(db[0])
        if (params.seq_type == 'pep') {
            blastp_nr(ch_split_seqs.nr, db, 'xml', '')
            parse_blastp_nr(blastp_nr.out.xml)
            parse_blastp_nr.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastp_nr_txt }
            combine_blastp_nr(blastp_nr_txt, 'blastp', entap_sequence_filename, 'nr')
        }
        if (params.seq_type == 'nuc') {
            blastx_nr(ch_split_seqs.nr, db, 'xml', '')
            parse_blastx_nr(blastx_nr.out.xml)
            parse_blastx_nr.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastx_nr_txt }
            combine_blastx_nr(blastx_nr_txt, 'blastx', entap_sequence_filename, 'nr')
        }
    }

    //
    // BLAST sequences against NCBI refseq using Diamond.
    //
    if (params.data_refseq) {
        db = [ file(params.data_refseq + "/refseq_plant.protein.dmnd", checkIfExists: true) ]
        entap_dbs.add(db[0])
        if (params.seq_type == 'pep') {
            blastp_refseq(ch_split_seqs.refseq, db, 'xml', '')
            parse_blastp_refseq(blastp_refseq.out.xml)
            parse_blastp_refseq.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastp_refseq_txt }
            combine_blastp_refseq(blastp_refseq_txt, 'blastp', entap_sequence_filename, 'refseq_plant')
        }
        if (params.seq_type == 'nuc') {
            blastx_refseq(ch_split_seqs.refseq, db, 'xml', '')
            parse_blastx_refseq(blastx_refseq.out.xml)
            parse_blastx_refseq.out.blast_txt
                .map { it[1] }
                .collect()
                .set { blastx_refseq_txt }
            combine_blastx_refseq(blastx_refseq_txt, 'blastx', entap_sequence_filename, 'refseq_plant')
        }
    }

    //
    // BLAST sequences against OrthoDB using Diamond.
    //
    if (params.data_orthodb) {
        db = [ file(params.data_orthodb, checkIfExists: true) ]
        if (params.seq_type == 'pep') {
            blastp_orthodb(ch_split_seqs.orthodb, db, 'xml', '')
            blastp_orthodb.out.txt
                .map { it[1] }
                .set { blastp_orthodb_txt }
            find_ortho_groups(blastp_orthodb_txt, sequence_filename)
        }
        if (params.seq_type == 'nuc') {
            blastx_orthodb(ch_split_seqs.orthodb, db, 'xml', '')
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
            blastp_string(ch_split_seqs.string, db, 'xml', '')
        }
        if (params.seq_type == 'nuc') {
            blastx_string(ch_split_seqs.string, db, 'xml', '')
        }
    }
    //
    // InterProScan
    //
    if (params.data_ipr) {

        if (params.seq_type == 'pep') {
            interproscan_pep(ch_split_seqs.ipr, params.seq_type)
            interproscan_pep.out.outfiles
                .map { it[1] }
                .flatten()
                .branch {
                    tsv: it.getFileName().toString().endsWith(".tsv")
                    xml: it.getFileName().toString().endsWith(".xml")
                }
                .set { interproscan_pep_out }
            interproscan_combine(interproscan_pep_out.tsv.collect(), sequence_filename)
        }
        if (params.seq_type == 'nuc') {
            interproscan_nuc(ch_split_seqs.ipr, params.seq_type)
            interproscan_nuc.out.outfiles
                .map { it[1] }
                .flatten()
                .branch {
                    tsv: it.getFileName().toString().endsWith(".tsv")
                    xml: it.getFileName().toString().endsWith(".xml")
                }
                .set { interproscan_nuc_out }
            interproscan_combine(interproscan_nuc_out.tsv.collect(), sequence_filename)
        }
    }

    //
    // EnTAP
    //
    if (!params.skip_entap) {
        entap_config("${baseDir}/assets/entap_config.ini")

        // Create a channel containing all of the combined blast output files
        ch_blast_results = Channel.empty()
        if (params.seq_type == 'pep') {
            if (params.data_sprot) {
                ch_blast_results = ch_blast_results.concat(combine_blastp_sprot.out.outfile);
            }
            if (params.data_refseq) {
                ch_blast_results = ch_blast_results.concat(combine_blastp_refseq.out.outfile);
            }
            if (params.data_nr) {
                ch_blast_results = ch_blast_results.concat(combine_blastp_nr.out.outfile);
            }
            if (params.data_trembl) {
                ch_blast_results = ch_blast_results.concat(combine_blastp_trembl.out.outfile);
            }
        }
        if (params.seq_type == 'nuc') {
            if (params.data_sprot) {
                ch_blast_results = ch_blast_results.concat(combine_blastx_sprot.out.outfile);
            }
            if (params.data_refseq) {
                ch_blast_results = ch_blast_results.concat(combine_blastx_refseq.out.outfile);
            }
            if (params.data_nr) {
                ch_blast_results = ch_blast_results.concat(combine_blastx_nr.out.outfile);
            }
            if (params.data_trembl) {
                ch_blast_results = ch_blast_results.concat(combine_blastx_trembl.out.outfile);
            }
        }

        entap_run(params.input, entap_dbs, "${baseDir}/assets/entap_config.ini",
            params.seq_type, entap_config.out.entap_db, entap_config.out.eggnog_db,
            entap_config.out.data_eggnog, ch_blast_results.collect(), interproscan_combine.out.tsv)
    }

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.adaptivecard(workflow, params, summary_params, projectDir, log)
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
