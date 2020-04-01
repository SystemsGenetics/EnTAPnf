#!/usr/bin/env python3

import os
import sys
import pandas as pd
from parse_blastxml import parseBLASTXMLfile
import sqlite3
from sqlite3 import Error

def create_connection(db_file):
    """
    Create a database connection to a SQLite database
    """

    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print("SQLite version: " + sqlite3.version)
    except Error as e:
        print("SQLite error: ")
        print(e)

    return conn


if __name__ == "__main__":
    """
    The main subroutine of the program.  The script receives a single arugment:
    the file name of the BLAST results against the OrthoDB FASTA file of
    protein sequences.
    """
    script_name, blast_results, orthodb_dir, outfile = sys.argv

    # Open a connection to the orthodb SQlite database.
    if not os.path.exists(orthodb_dir + "/odb10v1.db"):
        print("ERROR: The odb10v1.db index file cannot be found at {}.".format(orthodb_dir + "/odb10v1.db"))
        exit(1)

    # Connect to the SQLite database.
    conn = create_connection(orthodb_dir + "/odb10v1.db")

    # First parse the blast results.
    results = parseBLASTXMLfile(blast_results, True)

    # Reduce the list down to the best hit per query sequence.
    grp_q = results[['Query_id', 'Hit_id', 'Evalue', 'Identity']].groupby(by='Query_id')
    best_hits = grp_q.apply(lambda x: x[x['Evalue'] == x['Evalue'].max()].iloc[0])

    # Create the output dataframe
    outdf = pd.DataFrame(columns=['Gene_ID','OrthoDB_OG_ID','OG_Name','Tax_ID','Database', 'Term_ID', 'OG_Term_Genes'])

    for index,row in best_hits.iterrows():

        # Get the list of orthologous groups per match.
        sql = """
          SELECT OG_id FROM odb10v1_OG2genes WHERE odb_gene_id = ?;
        """
        cur = conn.cursor()
        cur.execute(sql, (row['Hit_id'],))
        OGs = cur.fetchall()

        # Iterate through the orthologous groups to get the inforamtion
        # such as the taxonomy ID, name and annotations.
        for OG in OGs:

            # First get the name and taxonomy ID for the group.
            sql = """
              SELECT ncbi_taxid, name FROM odb10v1_OGs WHERE OG_id = ?;
            """
            cur = conn.cursor()
            cur.execute(sql, (OG[0],))
            OG_info = cur.fetchall()[0]

            # Second get the annotations
            sql = """
              SELECT db_name, term_id, num_genes FROM odb10v1_OG_xrefs WHERE OG_id = ?;
            """
            cur = conn.cursor()
            cur.execute(sql, (OG[0],))
            annots = cur.fetchall()

            for annot in annots:
                outdf = outdf.append(pd.Series({
                  'Gene_ID': row['Query_id'],
                  'OrthoDB_OG_ID': OG[0],
                  'OG_Name': OG_info[1],
                  'Tax_ID': OG_info[0],
                  'Database': annot[0],
                  'Term_ID': annot[1],
                  'OG_Term_Genes': annot[2],
                }), ignore_index=True)

    # Finally, write the output file with the results.
    outdf.to_csv(outfile, sep="\t", index=False)
    exit(0)
