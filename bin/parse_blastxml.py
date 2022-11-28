#!/usr/bin/env python3

"""
Parses BLAST compatible XML files

.. module:: Annotater
    :platform: UNIX, Linux
    :synopsis: add synopsis.

"""

from lxml import etree
import pandas as pd
import argparse
import re


class HitTarget(object):
    """
    This class is used as a Target parser for the lxml XMLParser.
    """

    def __init__(self, parse_ncbi):
        """ "
        Constructor that sets the member variables.
        """
        self.parse_ncbi = parse_ncbi
        # An array of all the HSPs. This will get converted to a Pandas
        # data frame when all done.
        self.hit_list = []
        # This is set to True when we are in the Hit element.
        self.in_hit = False
        # This is set to True when we are in teh Hsp element.
        self.in_hsp = False
        # Holds a Series object for the current HSP. This gets refreshed
        # for every HSP.
        self.current_hsp = self.initHSP()
        # Stores the current tag that is being processed by the parser.
        self.current_tag = ""
        # These store details about the query sequence.
        self.current_query = ""
        self.current_query_def = ""
        self.current_query_len = ""
        # These store details about the hit sequence.
        self.current_hit = ""
        self.current_hit_def = ""
        self.current_hit_len = ""
        self.current_hit_acc = ""

    def initHSP(self):
        """ "
        Initializes a new Series object for a new HSP
        """
        return pd.Series(
            index=[
                "Query_id",
                "Query_def",
                "Query_len",
                "Hit_id",
                "Hit_def",
                "Hit_len",
                "Hit_accession",
                "Hsp_bit-score",
                "Hsp_score",
                "Hsp_evalue",
                "Hsp_identity",
                "Hsp_positive",
                "Hsp_gaps",
                "Hsp_align-len",
                "Hsp_query-from",
                "Hsp_query-to",
                "Hsp_hit-from",
                "Hsp_hit-to",
                "Hsp_query-frame",
                "Hsp_hit-frame",
            ],
            dtype=object,
        )

    def start(self, tag, attrib):
        """ "
        When a new tag is encountered in the XML this function is called.

        We can then use the tag to set state.

        :param tag: the XML tag that was encountered
        :param attrib: the list of attributes for the tag.
        """
        if tag == "Hit":
            self.in_hit = True
        if tag == "Hsp":
            self.in_hsp = True
        self.current_tag = tag

    def end(self, tag):
        """ "
        When a tag has been completely parsed this function is called.

        We can then use this to set state.

        :param tag: the XML tag that was encountered
        """
        if tag == "Hit":
            self.in_hit = False
        if tag == "Hsp":
            self.in_hsp = False
            # Add in the query details to the HSP
            self.current_hsp.loc["Query_id"] = self.current_query
            self.current_hsp.loc["Query_def"] = self.current_query_def
            self.current_hsp.loc["Query_len"] = self.current_query_len
            # Add in the hit details to the HSP
            if self.parse_ncbi:
                self.current_hsp.loc["Hit_id"] = self.current_hit
                self.current_hsp.loc["Hit_def"] = self.current_hit_def
                self.current_hsp.loc["Hit_len"] = self.current_hit_len
                self.current_hsp.loc["Hit_accession"] = self.current_hit_acc
            else:
                # TODO: adjust these for non NCBI recognized databases.

                self.current_hsp.loc["Hit_id"] = self.current_hit
                self.current_hsp.loc["Hit_def"] = self.current_hit_def
                self.current_hsp.loc["Hit_len"] = self.current_hit_len
                self.current_hsp.loc["Hit_accession"] = self.current_hit_acc
            # Add this HSP to the list.
            self.hit_list.append(self.current_hsp)
            self.current_hsp = self.initHSP()

    def data(self, data):
        """ "
        When a tag parsed this function is called to provide the data.

        If we are in an HSP tag we use this function to set the HSP attributes.

        :param data: the text data in the tag.
        """

        # Skip when the data is empty.
        if not data.strip():
            return

        # Handle query details
        if self.current_tag == "Iteration_query-ID":
            self.current_query = data

        if self.current_tag == "Iteration_query-def":
            self.current_query_def = data
            if self.parse_ncbi is None:
                match = re.match(r"^(.*?)\s.*$", self.current_query_def)
                if match:
                    self.current_query = match.group(1)

        if self.current_tag == "Iteration_query-len":
            self.current_query_len = data

        # Handle hit details
        if self.current_tag == "Hit_id":
            self.current_hit = data
        if self.current_tag == "Hit_def":
            self.current_hit_def = data
        if self.current_tag == "Hit_accession":
            self.current_hit_acc = data
        if self.current_tag == "Hit_len":
            self.current_hit_len = data

        # Ignore some fields but save all others if we're in an hsp
        ignore = ["Hsp_num", "Hsp_qseq", "Hsp_hseq", "Hsp_midline"]
        if self.current_tag in ignore:
            pass
        elif self.in_hsp:
            self.current_hsp.loc[self.current_tag] = data

    def close(self):
        """ "
        When the XML file is fully parsed this function is called.

        We use this function to build a final Pandas data frame containing
        all of the HSPs and to return it.
        """
        results = pd.DataFrame(self.hit_list)
        results = results.rename(
            columns={
                "Hsp_bit-score": "Bit_Score",
                "Hsp_score": "Score",
                "Hsp_evalue": "Evalue",
                "Hsp_identity": "Identity",
                "Hsp_positive": "Positive",
                "Hsp_gaps": "Gaps",
                "Hsp_align-len": "Align_len",
                "Hsp_query-from": "Query_from",
                "Hsp_query-to": "Query_to",
                "Hsp_hit-from": "Hit_from",
                "Hsp_hit-to": "Hit_to",
                "Hsp_query-frame": "Query_frame",
                "Hsp_hit-frame": "Hit_frame",
            }
        )

        # Remove some unwanted columns
        # results = results.drop('Query_def', 1)

        # Return the final resulting data frame.
        return results


def parseBLASTXMLfile(xml_file, parse_ncbi):
    """ "
    Parses a BLAST XML file and returns a list of hits in a Pandas Dataframe.

    :param xml_file:  The path (full or relative) to the XML file.
    """
    parser = etree.XMLParser(target=HitTarget(parse_ncbi))
    return etree.parse(xml_file, parser)


def main():
    """
    The main function
    """

    # Specifies the arguments for this script
    parser = argparse.ArgumentParser()
    parser.add_argument("--xml_file", dest="xml_file", type=str, required=True)
    parser.add_argument("--out_file", dest="out_file", type=str, required=True)
    parser.add_argument("--parse_ncbi", dest="parse_ncbi", type=bool, required=False)

    # Read in the input arguments
    args = parser.parse_args()

    # parse the XML file.
    results = parseBLASTXMLfile(args.xml_file, args.parse_ncbi)
    results.to_csv(args.out_file, sep="\t", index=False)


if __name__ == "__main__":
    main()
