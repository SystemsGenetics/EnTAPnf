#!/usr/bin/env python3
import sys
import argparse
import numpy as np
import pandas as pd
from parse_blastxml import parseBLASTXMLfile
from _enzymelookup import EnzymeLookup




parser = argparse.ArgumentParser()
parser.add_argument('--xml', dest='xmlPath', type=str, required=True)
parser.add_argument('--enzyme', dest='enzymePath', type=str, required=True)
parser.add_argument('--out', dest='outPath', type=str, required=True)
args = parser.parse_args()
data = parseBLASTXMLfile(args.xmlPath,True)
data = data[data['Evalue'].astype('float') <= 1e-50]
query_def = data["Query_def"]
hit_id = data["Hit_id"]
enzymeLookup = EnzymeLookup(args.enzymePath)

columns = ['Gene_ID','SwissProt_Accession','SwissProt_Entry_Name','EC_Number','EC_Description']
outdf = pd.DataFrame(columns=columns)

for (query,hit) in zip(query_def,hit_id):
    hit_parts = hit.split("|")
    rows = enzymeLookup.find(query.split()[0], hit_parts[1], hit_parts[2])
    for row in rows:
        outdf = outdf.append(pd.Series({
          'Gene_ID': row[0],
          'SwissProt_Accession': row[1],
          'SwissProt_Entry_Name': row[2],
          'EC_Number': row[3],
          'EC_Description': row[4],
        }), ignore_index=True)

outdf.drop_duplicates(inplace=True)
outdf.to_csv(args.outPath, sep="\t", index=False)
exit(0)
