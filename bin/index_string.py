#!/usr/bin/env python3
"""
This is the index string script.
"""
import argparse
import sqlite3


#
# The expected number of columns for each row of the protein links CSV file.
#
LINKS_CSV_COL_SIZE = 16


#
# The expected number of columns for each row of the protein info CSV file.
#
INFO_CSV_COL_SIZE = 4


#
# The extension used for the protein info database.
#
INFO_DB_EXT = ".info.db"


#
# The extension used for all the protein links databases.
#
LINKS_DB_EXT = ".links.db"


#
# The maximum number of rows inserted into one of the multiple databases for
# protein links data.
#
ROWS_PER_DB = 620000000




def makeInfoDB(
    path
    ):
    """
    Creates a new database and the required table for protein info data.

    Parameters
    ----------
    path : string
           Path of the new SQL database that is created.

    Returns
    -------
    ret0 : sqlite3.Connection
           A connection to the new database.
    """
    if os.path.exists(path):
        os.remove(path)
    conn = sqlite3.connect(path)
    sql = (
        "CREATE TABLE IF NOT EXISTS protein_info ("
        "protein_external_id text NOT NULL,"
        "preferred_name text NOT NULL,"
        "protein_size integer NOT NULL,"
        "annotation text NOT NULL);"
    )
    conn.cursor().execute(sql)
    return conn




def makeLinksDB(
    path
    ):
    """
    Creates a new database and the required table for protein links data.

    Parameters
    ----------
    path : string
           Path of the new SQL database that is created.

    Returns
    -------
    ret0 : sqlite3.Connection
           A connection to the new database.
    """
    conn = sqlite3.connect(path)
    sql = (
        "CREATE TABLE IF NOT EXISTS protein_links_full ("
        "protein1 text NOT NULL,"
        "protein2 text NOT NULL,"
        "neighborhood integer NOT NULL,"
        "neighborhood_transferred integer NOT NULL,"
        "fusion integer NOT NULL,"
        "cooccurence integer NOT NULL,"
        "homology integer NOT NULL,"
        "coexpression integer NOT NULL,"
        "coexpression_transferred integer NOT NULL,"
        "experiments integer NOT NULL,"
        "experiments_transferred integer NOT NULL,"
        "database integer NOT NULL,"
        "database_transferred integer NOT NULL,"
        "textmining integer NOT NULL,"
        "textmining_transferred integer NOT NULL,"
        "combined_score integer NOT NULL);"
    )
    conn.cursor().execute(sql)
    return conn




def populateInfoDB(
    conn
    ,ifile
    ):
    """
    Populates the given database with protein info data from the entire contents
    of the given CSV input file.

    Parameters
    ----------
    conn : sqlite3.Connection
           A connection to the SQL database that is populated with protein info
           data.
    ifile : io.TextIOWrapper
            The text CSV protein info file where data is pulled from to populate
            the database.
    """
    while True:
        line = ifile.readline()
        if not line:
            conn.cursor().execute(
                "CREATE INDEX idx_protein_info_protein_external_id"
                " ON protein_info (protein_external_id);"
            )
            conn.commit()
            break
        parts = line[:-1].split("\t")
        if len(parts)==INFO_CSV_COL_SIZE:
            conn.cursor().execute(
                "INSERT INTO protein_info"
                "(protein_external_id,preferred_name,protein_size,annotation)"
                "VALUES (?,?,?,?);"
                ,parts
            )




def populateLinksDB(
    conn
    ,ifile
    ):
    """
    Populates the given database with protein links data from the given CSV
    input file. After a certain number of rows have been inserted or the end of
    the CSV input file has been reached it created appropriate indexes and
    returns.

    Parameters
    ----------
    conn : sqlite3.Connection
           A connection to the SQL database that is populated with protein links
           data.
    ifile : io.TextIOWrapper
            The text CSV protein links full file where data is pulled from to
            populate the database.

    Returns
    -------
    ret0 : bool
           True if there is more data to read from the given input file or false
           if the end of file has been reached and there is no more data to
           populate to the next database.
    """
    sqlIdx1 = (
        "CREATE INDEX idx_protein_links_full_protein1"
        " ON protein_links_full (protein1);"
    )
    sqlIdx2 = (
        "CREATE INDEX idx_protein_links_full_protein2"
        " ON protein_links_full (protein2);"
    )
    count = 0
    while True:
        line = ifile.readline()
        if not line:
            conn.cursor().execute(sqlIdx1)
            conn.cursor().execute(sqlIdx2)
            conn.commit()
            return False
        parts = line.split()
        if len(parts)==LINKS_CSV_COL_SIZE:
            conn.cursor().execute(
                "INSERT INTO protein_links_full"
                "(protein1,protein2,neighborhood,neighborhood_transferred,fusion,cooccurence"
                ",homology,coexpression,coexpression_transferred,experiments"
                ",experiments_transferred,database,database_transferred,textmining"
                ",textmining_transferred,combined_score)"
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?);"
                ,parts
            )
            count += 1
            if count==ROWS_PER_DB:
                conn.cursor().execute(sqlIdx1)
                conn.cursor().execute(sqlIdx2)
                conn.commit()
                return True








### __SCRIPT__ ###
parser = argparse.ArgumentParser()
parser.add_argument('--links',dest='linksPath',type=str,required=True)
parser.add_argument('--info',dest='infoPath',type=str,required=True)
parser.add_argument('--out',dest='outPath',type=str,required=True)
args = parser.parse_args()
conn = makeInfoDB(args.outPath+INFO_DB_EXT)
with open(args.infoPath,"r") as ifile:
    ifile.readline()
    populateInfoDB(conn,ifile)
dbIdx = 0
conn = makeLinksDB("%s%02d%s"%(args.outPath,dbIdx,LINKS_DB_EXT))
with open(args.linksPath,"r") as ifile:
    ifile.readline()
    while populateLinksDB(conn,ifile):
        dbIdx += 1
        conn = makeLinksDB("%s%02d%s"%(args.outPath,dbIdx,LINKS_DB_EXT))
