//
// This file holds several functions specific to the workflow/entap.nf in the nf-core/entap pipeline
//
import java.io.File

class WorkflowEntap {

    //
    // Check and validate parameters
    //
    public static void initialise(params, log) {
        File data_file
        File data_dir

        // Make sure that the SwissProt BLAST settings are good.
        if (params.data_sprot) {
          // Make sure the data directory is present.
          data_file = new File("${params.data_sprot}/uniprot_sprot.dmnd")
          if (data_file.isEmpty()) {
            log.error "Error: the Uniprot SwissProt Diamond index file cannot be found at ${params.data_sprot}/uniprot_sprot.dmnd. Please check the params.data_sprot setting and make sure the file is present in the specified directory."
          }
        }

        // Make sure that if the NR BLAST settings are good.
        if (params.data_nr) {
          // Make sure the data directory is present.
          data_file = new File("${params.data_nr}/nr.dmnd")
          if (data_file.isEmpty()) {
            log.error "Error: the NCBI nr Diamond index file cannot be found at ${params.data_nr}/nr.dmnd. Please check the params.data_nr setting and make sure the file is present in the specified directory."
          }
        }

        // Make sure that if the NR BLAST settings are good.
        if (params.data_refseq) {
          // Make sure the data directory is present.
          data_file = new File("${params.data_refseq}/refseq_plant.protein.dmnd")
          if (data_file.isEmpty()) {
            log.error "Error: the NCBI nr RefSeq index file cannot be found at ${params.data_refseq}/refseq_plant.protein.dmnd. Please check the params.data_refseq setting and make sure the file is present in the specified directory."
          }
        }

        // Make sure that if the Trembl BLAST settings are good.
        if (params.data_trembl) {
          // Make sure the data directory is present.
          data_file = new File("${params.data_trembl}/uniprot_trembl.dmnd")
          if (data_file.isEmpty()) {
            log.error "Error: the Uniprot trembl Diamond index file cannot be found at ${params.data_trembl}/uniprot_trembl.dmnd. Please check the params.data_trembl setting and make sure the file is present in the specified directory."
          }
        }

        // Make sure that the String BLAST settings are good.
        if (params.data_string) {
          // Make sure the data directory is present.
          data_file = new File("${params.data_string}/protein.sequences.v11.0.dmnd")
          if (data_file.isEmpty()) {
            log.error "Error: the String Diamond index file cannot be found at ${params.data_string}/protein.sequences.v11.0.dmnd. Please check the params.data_string setting and make sure the file is present in the specified directory."
          }
          data_file = new File("${params.data_string}/protein.links.full.v11.0.txt")
          if (data_file.isEmpty()) {
            log.error "Error: the String linksfile cannot be found at ${params.data_string}/protein.links.full.v11.0.txt. Please check the params.data_string setting and make sure the file is present in the specified directory."
          }
          data_file = new File("${params.data_string}/protein.info.v11.0.txt")
          if (data_file.isEmpty()) {
            log.error "Error: the String Diamond index file cannot be found at ${params.data_string}/protein.info.v11.0.txt. Please check the params.data_string setting and make sure the file is present in the specified directory."
          }
        }

        // Make sure that if the OrthoDB database is specified the settings are good.
        if (params.data_orthodb) {
          // Make sure the data directory is present.
          data_dir = new File("${params.data_orthodb}")
          if (data_dir.isEmpty()) {
            log.error "Error: the OrthoDB data directory cannot be found: ${params.data_orthodb}. Please check the params.data_orthodb setting."
          }
        }

        // Make sure that if the InterProScan settings are good.
        if (params.data_ipr) {
          // TODO: make sure the applicatin list is valid.

          // Make sure the data directory is present.
          data_dir = new File("${params.data_ipr}")
          if (data_dir.isEmpty()) {
            log.error "Error: the InterProScan data directory cannot be found: ${params.data_ipr}. Please check the params.data_ipr setting."
          }
        }
    }
}
