#!/usr/bin/env python3
"""
Retrieves all species below a given taxonomic level in OrthoDB.

.. module:: Annotater
    :platform: UNIX, Linux
    :synopsis: add synopsis.

"""
import argparse
import sys
import pandas as pd
import numpy as np
import re


def main():
    # Specifies the arguments for this script
    parser = argparse.ArgumentParser()
    parser.add_argument('taxid', action='store')
    parser.add_argument('level2species', action='store')

    # Read in the input arguments
    args = parser.parse_args()

    # Maps all of the species under each level of the OrthoDB tree.
    mapping = {}

    # Open the level2species file and store all of the species for each
    # taxonomic level for easy lookup
    col_names = np.array(['Top_Level', 'Org_ID', 'Num_Hops', 'Levels'])
    level2species = pd.read_csv(args.level2species, sep='\t', header=None, names=col_names)
    for index, row in level2species.iterrows():
        if (row['Num_Hops'] > 1):

            # Remove the _0 from the taxonomic ID and parens from the levesl
            levels = row['Levels']
            org_id = row['Org_ID']
            levels = levels.replace('}', '')
            levels = levels.replace('{', '')
            levels = levels.split(',')

            # Iterate through the levels and add the species on the end of the
            # list of levels to every parent in the list.
            for i in range(0, len(levels)-1):
                parent = levels[i]
                species = levels[len(levels)-1]
                if parent not in mapping:
                    mapping[parent] = {}
                if species not in mapping:
                    mapping[species] = {}
                mapping[parent][species] = org_id
                mapping[species][species] = org_id

    # Now print out the keys
    print(",".join(mapping[args.taxid].values()))

if __name__ == "__main__":
    main()
