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


#
# The root path where all database files for indexed protein info and protein
# links can be found. The appropriate extension is added as required for opening
# any specific database. This is set in the script section from a command line
# argument.
#
databasePath = ""




def filterBlastXML(
    path
    ):
    """
    Getter function.

    Parameters
    ----------
    path : object
           Path to the blast XML file that is filtered.

    Returns
    -------
    ret0 : list
           Filtered proteins from the blast XML file at the given path. The
           proteins are filtered based off finding each unique protein with the
           highest correlation length or e-value.
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
    protein
    ):
    """
    Getter method.

    Parameters
    ----------
    protein : object
              Protein id that is matched with protein link rows protein columns.

    Returns
    -------
    ret0 : list
           Protein link database rows that match the given protein with the rows
           protein1 or protein2 column. This queries all separate protein link
           databases.
    """
    sql = "SELECT * FROM protein_links_full WHERE `protein1` = ? OR `protein2` = ?;"
    i = 0
    ret = []
    while True:
        rpath = "%s%02d.links.db"%(databasePath,i)
        if os.path.isfile(rpath):
            conn = sqlite3.connect(rpath)
            csr = conn.cursor()
            csr.execute(sql,(protein,protein))
            ret += csr.fetchall()
        else:
            return ret
        i += 1




def processProteinLinks(
    graph
    ,rows
    ):
    """
    Processes the given list of protein link database rows, adding their
    vertices and edges to the given igraph object. All required attributes for
    proteins and their links are added to their respective vertices and edges in
    the igraph object. The order of the two protein ids for edges does not
    matter so reversed links are ignored as duplicates.

    Parameters
    ----------
    graph : igraph.Graph
            The graph object that is populated with protein vertices and protein
            link edges.
    rows : list
           The protein link database rows that is used to populate the given
           igraph object.
    """
    proteins = set()
    links = set()
    def addVertex(protein):
        if protein not in proteins:
            info = proteinInfo(protein)
            if info is not None:
                graph.add_vertices(protein)
                vertex = graph.vs[-1]
                vertex["preferred_name"] = info[1]
                vertex["protein_size"] = info[2]
                vertex["annotation"] = info[3]
                proteins.add(protein)
    for row in rows:
        addVertex(row[0])
        addVertex(row[1])
        edgeName = row[0]+"__"+row[1]
        altEdgeName = row[1]+"__"+row[0]
        if edgeName not in links and altEdgeName not in links:
            graph.add_edges(((row[0],row[1]),))
            edge = graph.es[-1]
            edge["name"] = edgeName
            edge["neighborhood"] = row[2]
            edge["neighborhood_transferred"] = row[3]
            edge["fusion"] = row[4]
            edge["cooccurence"] = row[5]
            edge["homology"] = row[6]
            edge["coexpression"] = row[7]
            edge["coexpression_transferred"] = row[8]
            edge["experiments"] = row[9]
            edge["experiments_transferred"] = row[10]
            edge["database"] = row[11]
            edge["database_transferred"] = row[12]
            edge["textmining"] = row[13]
            edge["textmining_transferred"] = row[14]
            edge["combined_score"] = row[15]
            links.add(edgeName)




def proteinInfo(
    protein
    ):
    """
    Getter method.

    Parameters
    ----------
    protein : object
              The external protein id that is used to find a protein info
              database match.

    Returns
    -------
    ret0 : object
           A single protein info database row that matches the given protein id
           or None if no match is found in the protein info database.
    """
    sql = "SELECT * FROM protein_info WHERE `protein_external_id` = ?;"
    rpath = "%s.info.db"%(databasePath,)
    conn = sqlite3.connect(rpath)
    csr = conn.cursor()
    csr.execute(sql,(protein,))
    ret = csr.fetchall()
    return ret[0] if ret else None








### __SCRIPT__ ###
parser = argparse.ArgumentParser()
parser.add_argument('--xml',dest='xmlPath',type=str,required=True)
parser.add_argument('--db',dest='dbPath',type=str,required=True)
parser.add_argument('--species',dest='speciesId',type=str,required=True)
parser.add_argument('--out',dest='outPath',type=str,required=True)
args = parser.parse_args()
databasePath = args.dbPath
results = filterBlastXML(args.xmlPath)
graph = igraph.Graph()
rows = []
for protein in results:
    rows += findMatches(protein)
processProteinLinks(graph,rows)
with open(args.outPath,"w") as ofile:
    graph.write(ofile,"graphml")
