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

class HitTarget(object):
    """
    This class is used as a Target parser for the lxml XMLParser.
    """

    def __init__(self):
        """"
        Constructor that sets the member variables.
        """
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
        self.current_tag = ''

    def initHSP(self):
        """"
        Initializes a new Series object for a new HSP
        """
        return pd.Series(index=['Hit_num', 'Hsp_num', 'Hsp_bit-score', 'Hsp_score', 'Hsp_evalue', 'Hsp_query-from', 'Hsp_query-to', 'Hsp_hit-from', 'Hsp_hit-to', 'Hsp_query-frame', 'Hsp_identity', 'Hsp_positive', 'Hsp_gaps', 'Hsp_align-len', 'Hsp_hseq', 'Hsp_midline'])

    def start(self, tag, attrib):
        """"
        When a new tag is encountered in the XML this function is called.

        We can then use the tag to set state.

        :param tag: the XML tag that was encountered
        :param attrib: the list of attributes for the tag.
        """
        if tag == 'Hit':
            self.in_hit = True
        if tag == 'Hsp':
            self.in_hsp = True
        self.current_tag = tag

    def end(self, tag):
        """"
        When a tag has been completely parsed this function is called.

        We can then use this to set state.

        :param tag: the XML tag that was encountered
        """
        if tag == 'Hit':
            self.in_hit = False
        if tag == 'Hsp':
            self.in_hsp = False
            self.hit_list.append(self.current_hsp)
            self.current_hsp = self.initHSP()


    def data(self, data):
        """"
        When a tag parsed this function is called to provide the data.

        If we are in an HSP tag we use this function to set the HSP attributes.

        :param data: the text data in the tag.
        """

        # Skip when the data is empty.
        if (not data.strip()):
            return;

        # Handle the HSP tag names.
        if self.current_tag == 'Hit_num':
            self.current_hsp.loc['Hit_num'] = data
        elif self.in_hsp:
            self.current_hsp.loc[self.current_tag] = data


    def close(self):
        """"
        When the XML file is fully parsed this function is called.

        We use this function to build a final Pandas data frame containing
        all of the HSPs and to return it.
        """
        return(pd.DataFrame(self.hit_list))

def parseBLASTXMLfile(xml_file):
    """"
    Parses a BLAST XML file and returns a list of hits in a Pandas Dataframe.

    :param xml_file:  The path (full or relative) to the XML file.
    """
    parser = etree.XMLParser(target = HitTarget())
    return etree.parse(xml_file, parser)


def main():
    """
    The main function
    """

    # Specifies the arguments for this script
    parser = argparse.ArgumentParser()
    parser.add_argument('xml_file', action='store')

    # Read in the input arguments
    args = parser.parse_args()

    # parse the XML file.
    print(parseBLASTXMLfile(args.xml_file))

if __name__ == "__main__":
    main()
