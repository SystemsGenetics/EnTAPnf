#!/usr/bin/env python3
"""
This is the parse enzyme script.
"""
import sys
import argparse
import numpy as np
import pandas as pd
from parse_blastxml import parseBLASTXMLfile
from _enzymelookup import EnzymeLookup








class EnzymeLookup():
    """
    This is the enzyme lookup class. It is a reverse lookup table that loads an
    enzyme data from the path given to its initialize operator. In turn a find
    method is supplied that looks up any matching names.
    """


    #######################
    # PUBLIC - Initialize #
    #######################


    def __init__(
        self
        ,path
        ):
        """
        Initializes this new enzyme lookup table with the given path to the
        enzyme data file. The enzyme data file's ID, DE, and DR lines are loaded
        into this new table.

        Parameters
        ----------
        path : string
               File path to the enzyme data file this new table loads into
               memory.
        """
        self.__links = {}
        self.__loadReverse_(path)


    ####################
    # PUBLIC - Methods #
    ####################


    def find(
        self
        ,gene
        ,id_
        ,name
        ):
        """
        Getter method.

        Parameters
        ----------
        gene : string
               The gene that is matched. This is only used to add to returned
               tuples and not matched against.
        id_ : string
              The protein id that is matched. This is only used to add to
              returned tuples and not matched against.
        name : string
               The protein name that is matched. This is the value actually used
               to reverse lookup any EC numbers and descriptions that contain
               this protein name.

        Returns
        -------
        ret0 : list
               Tuples that contain any matches found. Each tuple contains the
               given gene, given protein id, given protein name, matched EC id,
               and matched EC description.
        """
        ret = []
        if name in self.__links:
            for node in self.__links[name]:
                ret.append((gene, id_, name, node["id"], node["de"]))
        return ret


    #####################
    # PRIVATE - Methods #
    #####################


    def __loadReverse_(
        self
        ,path
        ):
        """
        Loads this instance's lookup table from the given enzyme data file path.

        Parameters
        ----------
        path : string
               File path to the enzyme data file this new table loads into
               memory.
        """
        with open(path,"r") as ifile:
            READ_NODE = 0
            READ_LINKS = 1
            state = READ_LINKS
            node = None
            while True:
                line = ifile.readline()
                if not line:
                    break
                if state == READ_NODE:
                    if line[:2] == "DE":
                        node["de"] = line[5:-1]
                        state = READ_LINKS
                elif state == READ_LINKS:
                    if line[:2] == "DR":
                        for segment in line[5:-1].split():
                            if segment.endswith(";"):
                                l = self.__links.get(segment[:-1],[])
                                l.append(node)
                                self.__links[segment[:-1]] = l
                    elif line[:2] == "ID":
                        node = {"id": line[5:-1], "de": None}
                        state = READ_NODE








### __SCRIPT__ ###
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
