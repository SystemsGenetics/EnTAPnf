#!/usr/bin/env python3
"""
This is the parse protein protein script.
"""
import argparse
import os
import sqlite3
import sys
import igraph
from parse_blastxml import parseBLASTXMLfile




def filterBlastXML(
    path
    ):
    """
    Detailed description.

    Parameters
    ----------
    path : object
           Detailed description.
    """
    data = parseBLASTXMLfile(path,True)
    queryDefs = data["Query_def"]
    hitIds = data["Hit_id"]
    hitLens = data["Hit_len"]
    eValues = data["Evalue"]
    ret = {}
    for (queryDef,hitId,hitLen,eValue) in zip(queryDefs,hitIds,hitLens,eValues):
        parts = hitId.split(".")
        if parts and parts[0]==args.speciesId:
            #hitId = ".".join(parts[1:])
            if queryDef not in ret:
                ret[queryDef] = (hitId,hitLen,eValue)
            elif (
                hitLen>ret[queryDef][1]
                or (hitLen==ret[queryDef][1] and eValue>ret[queryDef][2])
            ):
                ret[queryDef] = (hitId,hitLen,eValue)
            elif hitLen==ret[queryDef][1] and eValue==ret[queryDef][2]:
                raise RuntimeError(
                    "ERROR: Blast rows found with identical alignment length and e-value!"
                )
    return set((r[0] for r in ret.values()))




def findMatches(
    path
    ,protein
    ):
    """
    Detailed description.

    Parameters
    ----------
    path : object
           Detailed description.
    protein : object
              Detailed description.
    """
    sql = "SELECT * FROM protein_links_full WHERE `protein1` = ? OR `protein2` = ?;"
    i = 0
    ret = []
    while True:
        rpath = "%s%02d.links.db"%(path,i)
        if os.path.isfile(rpath):
            conn = sqlite3.connect(rpath)
            csr = conn.cursor()
            csr.execute(sql,(protein,protein))
            ret += csr.fetchall()
        else:
            return ret
        i += 1








### __SCRIPT__ ###
parser = argparse.ArgumentParser()
parser.add_argument('--xml',dest='xmlPath',type=str,required=True)
parser.add_argument('--db',dest='dbPath',type=str,required=True)
parser.add_argument('--species',dest='speciesId',type=str,required=True)
args = parser.parse_args()
results = filterBlastXML(args.xmlPath)
g = igraph.Graph()
for protein in results:
    print("\n\n===============================================================================")
    print(protein)
    print("")
    for row in findMatches(args.dbPath,protein):
        if not g:
            g.add_vertices(row[0])
            g.add_vertices(row[1])
            edgeName = row[0]+"<->"+row[1]
            g.add_edges(((row[0],row[1]),))
            g.es[:-1]["name"] = edgeName
        else:
            if not g or not g.vs.select(name=row[0]):
                g.add_vertices(row[0])
            if not g.vs.select(name=row[1]):
                g.add_vertices(row[1])
            edgeName = row[0]+"<->"+row[1]
            if not g.es.select(name=edgeName):
                g.add_edges(((row[0],row[1]),))
                g.es[:-1]["name"] = edgeName
        print(row)
print("\n\n")
g.write(sys.stdout,"graphml")
