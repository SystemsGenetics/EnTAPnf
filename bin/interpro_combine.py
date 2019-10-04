#!/usr/bin/env python3


"""
Overall description of what this module does.

.. module:: Annotater
    :platform: UNIX, Linux
    :synopsis: add a short synoposis.

"""

import glob
import pandas as pd
import numpy as np
import re

def write_IPR(ipr_file, tsv_file):
    """"
    Retrieves the IPR terms from a TSV file and writes them to a mapping file.

    :param ipr_file: the file handle to the IPR_Mapping.txt file.
    :param tsv_file: the TSV file from InterProScan.

    """
    # The IPR results do not have a consistent number of columns. Here we creaqte
    # an array that gets used in the read_table function below that forces
    # the function to include 15 columns.
    cols = np.arange(15)

    # Import partial IPR results for this file and set the column anames,
    # drop rows that has NaN in IPR column, and remove duplicates.
    ipr_terms = pd.read_csv(tsv_file, sep='\t', header=None, names=cols, usecols=[0,11,12])
    ipr_terms.columns = ["Gene", "IPR", "Description"]
    ipr_terms = ipr_terms.dropna(subset=["IPR"])
    ipr_terms = ipr_terms.drop_duplicates(keep='first')

    # Removes the training ID number that INterProScan adds to the sequence
    ipr_terms["Gene"] = ipr_terms["Gene"].str.replace(r'^(.*?)_\d+', r'\1')
    ipr_terms.to_csv(ipr_file, sep="\t", mode='a', header=False, index=False)


def write_GO(go_file, tsv_file):
    """"
    Retrieves the GO terms from a TSV file and writes them to a mapping file.

    :param go_file: the file handle to the GO_Mapping.txt file.
    :param tsv_file: the TSV file from InterProScan.
    """

    # The IPR results do not have a consistent number of columns. Here we creaqte
    # an array that gets used in the read_table function below that forces
    # the function to include 15 columns.
    cols = np.arange(15)

    # Import partial GO results for this file and set the column anames,
    # drop rows that has NaN in IPR column, and remove duplicates.
    go_terms = pd.read_csv(tsv_file,sep='\t',header=None,names=cols,usecols=[0,13])
    go_terms.columns = ["Gene", "GO"]
    go_terms = go_terms.dropna(subset=["GO"])
    go_terms = go_terms.drop_duplicates(keep='first')

    # Removes the training ID number that INterProScan adds to the sequence
    go_terms["Gene"] = go_terms["Gene"].str.replace(r'^(.*?)_\d+', r'\1')

    print(go_terms)


if __name__ == "__main__":

    # Get a list of all TSV files.
    tsv_filenames = sorted(glob.glob('*.tsv'))

    # Open the IPR_mappintgs.txt file where the final IPR mappings will be stored.
    ipr_file = open('IPR_mappings.txt','w')
    go_file = open('GO_mappings.txt','w')

    # Set the header names for the IPR_mappings.txt file
    ipr_headers = ["Gene", "IPR", "Description"]
    ipr_file.write("\t".join(ipr_headers))
    ipr_file.write("\n")

    # Set the header names for the GO_mappings.txt file
    go_headers = ["Gene", "GO"]
    go_file.write("\t".join(go_headers))
    go_file.write("\n")

    # Iterate through each of the TSV files and pull out the IPR mappings.
    for tsv_file in tsv_filenames:
       write_IPR(ipr_file, tsv_file)
       write_GO(go_file, tsv_file)

    ipr_file.close()
    go_file.close()
