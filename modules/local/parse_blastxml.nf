
process PARSE_BLASTXML {
    publishDir "${params.outdir}/interproscan_combined",
        mode: params.publish_dir_mode

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
