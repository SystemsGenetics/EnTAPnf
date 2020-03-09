import sys
import argparse
from parse_blastxml import parseBLASTXMLfile
from _enzymelookup import EnzymeLookup




parser = argparse.ArgumentParser()
parser.add_argument('--xml', dest='xmlPath', type=str, required=True)
parser.add_argument('--enzyme', dest='enzymePath', type=str, required=True)
parser.add_argument('--out', dest='outPath', type=str, required=True)
args = parser.parse_args()
data = parseBLASTXMLfile(args.xmlPath,True)
query_def = data["Query_def"]
hit_id = data["Hit_id"]
enzymeLookup = EnzymeLookup(args.enzymePath)
with open(args.outPath,"w") as ofile:
    for (query,hit) in zip(query_def,hit_id):
        hit_parts = hit.split("|")
        rows = enzymeLookup.find(query.split()[0],hit_parts[1],hit_parts[2])
        for row in rows:
            print("%s\t%s\t%s\t%s\t%s" % row,file=ofile)
