#!/usr/bin/env python3
"""
This is the parse protein protein script.
"""
import argparse
import os
import sqlite3
import sys
from parse_blastxml import parseBLASTXMLfile








### __SCRIPT__ ###
def find_matches(path,protein):
    sql = "SELECT * FROM protein_links_full WHERE `protein1` = ? OR `protein2` = ?;"
    i = 0
    while True:
        rpath = "%s%03d.db"%(path,i)
        if os.path.isfile(rpath):
            conn = sqlite3.connect(rpath)
            csr = conn.cursor()
            csr.execute(sql,(protein,protein))
            for row in csr.fetchall():
                print(row)
        else:
            break
        i += 1


parser = argparse.ArgumentParser()
parser.add_argument('--xml',dest='xmlPath',type=str,required=True)
parser.add_argument('--db',dest='dbPath',type=str,required=True)
parser.add_argument('--species',dest='speciesId',type=str,required=True)
args = parser.parse_args()


data = parseBLASTXMLfile(args.xmlPath,True)
queryDefs = data["Query_def"]
hitIds = data["Hit_id"]
hitLens = data["Hit_len"]
eValues = data["Evalue"]
results = {}
for (queryDef,hitId,hitLen,eValue) in zip(queryDefs,hitIds,hitLens,eValues):
    parts = hitId.split(".")
    if parts and parts[0]==args.speciesId:
        if queryDef not in results:
            results[queryDef] = (hitId,hitLen,eValue)
        elif (
            hitLen>results[queryDef][1]
            or (hitLen==results[queryDef][1] and eValue>results[queryDef][2])
        ):
            results[queryDef] = (hitId,hitLen,eValue)
        elif hitLen==results[queryDef][1] and eValue==results[queryDef][2]:
            print(
                "ERROR: Blast rows found with identical alignment length and e-value!"
                ,file=sys.stderr
            )
            sys.exit(1)


for queryDef in results:
    print("\n\n===============================================================================")
    print(queryDef,results[queryDef])
    print("")
    find_matches(args.dbPath,results[queryDef][0])
