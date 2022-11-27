#!/usr/bin/env python3
"""
A Python script for indexing the OrthoDB data files for quick lookup by
AnnoTater.

.. module:: Annotater
    :platform: UNIX, Linux
    :synopsis: add synopsis.

"""

import sqlite3
from sqlite3 import Error
import re
import os
import sys


def create_connection(db_file):
    """create a database connection to a SQLite database"""
    conn = None
    try:
        conn = sqlite3.connect(db_file)
        print("SQLite version: " + sqlite3.version)
    except Error as e:
        print(e)

    return conn


def create_table(conn, sql):
    """create a table from the create_table_sql statement
    :param conn: Connection object
    :param sql: a CREATE TABLE statement
    :return:
    """
    try:
        c = conn.cursor()
        c.execute(sql)
    except Error as e:
        print(e)


def index_odb10v1_OG_xrefs(conn, orthodb_dir):
    """Creates and inserts data into the odb10v1_OG_xrefs table.
    :param conn: Connection object
    :param orthodb_dir: The directory where the orthodb files are kept.
    :return:
    """

    print("Indexing odb10v1_OG_xrefs.tab")

    # Create the OG_xrefs table.
    sql = """
      CREATE TABLE IF NOT EXISTS odb10v1_OG_xrefs (
        OG_id integer NOT NULL,
        db_name text NOT NULL,
        term_id text NOT NULL,
        num_genes integer NOT NULL);
    """
    create_table(conn, sql)

    # Add the data.
    sql = """
      INSERT INTO odb10v1_OG_xrefs(OG_id, db_name, term_id, num_genes)
      VALUES(?,?,?,?);
    """
    fh = open(orthodb_dir + "/odb10v1_OG_xrefs.tab")
    for line in fh:
        cols = line.strip("\n").split("\t")
        cur = conn.cursor()
        cur.execute(sql, cols)

    # Create the indexes.
    sql = "CREATE INDEX odb10v1_OG_xrefs_idx1 ON odb10v1_OG_xrefs(OG_id);"
    cur.execute(sql)

    # Commit the transaction.
    conn.commit()


def index_odb10v1_OGs(conn, orthodb_dir):
    """Creates and inserts data into the odb10v1_OGs table.
    :param conn: Connection object
    :param orthodb_dir: The directory where the orthodb files are kept.
    :return:
    """

    print("Indexing odb10v1_OGs.tab")

    # Create the OGs table.
    sql = """
      CREATE TABLE IF NOT EXISTS odb10v1_OGs (
        OG_id integer NOT NULL,
        ncbi_taxid integer NOT NULL,
        name text NOT NULL);
    """
    create_table(conn, sql)

    # Add the data.
    sql = """
      INSERT INTO odb10v1_OGs(OG_id, ncbi_taxid, name)
      VALUES(?,?,?);
    """
    fh = open(orthodb_dir + "/odb10v1_OGs.tab")
    for line in fh:
        cols = line.strip("\n").split("\t")
        cur = conn.cursor()
        cur.execute(sql, cols)

    # Create the indexes.
    sql = "CREATE INDEX odb10v1_OGs_idx1 ON odb10v1_OGs(OG_id);"
    cur.execute(sql)
    sql = "CREATE INDEX odb10v1_OGs_idx2 ON odb10v1_OGs(ncbi_taxid);"
    cur.execute(sql)

    # Commit the transaction.
    conn.commit()


def index_odb10v1_OG2genes(conn, orthodb_dir):
    """Creates and inserts data into the odb10v1_OG2genes table.
    :param conn: Connection object
    :param orthodb_dir: The directory where the orthodb files are kept.
    :return:
    """

    print("Indexing odb10v1_OG2genes.tab")

    # Create the OG2genes table.
    sql = """
      CREATE TABLE IF NOT EXISTS odb10v1_OG2genes (
        OG_id integer NOT NULL,
        odb_gene_id integer NOT NULL);
    """
    create_table(conn, sql)

    # Add the data.
    sql = """
      INSERT INTO odb10v1_OG2genes(OG_id, odb_gene_id)
      VALUES(?,?);
    """
    fh = open(orthodb_dir + "/odb10v1_OG2genes.tab")
    for line in fh:
        cols = line.strip("\n").split("\t")
        cur = conn.cursor()
        cur.execute(sql, cols)

    # Create the indexes.
    sql = "CREATE INDEX odb10v1_OG2genes_idx1 ON odb10v1_OG2genes(OG_id);"
    cur.execute(sql)
    sql = "CREATE INDEX odb10v1_OG2genes_idx2 ON odb10v1_OG2genes(odb_gene_id);"
    cur.execute(sql)

    # Commit the transaction.
    conn.commit()


def index_odb10v1_level2species(conn, orthodb_dir):
    """Creates and inserts data into the odb10v1_level2species table.
    :param conn: Connection object
    :param orthodb_dir: The directory where the orthodb files are kept.
    :return:
    """

    print("Indexing odb10v1_level2species.tab")

    # Create the level2species table.
    sql = """
      CREATE TABLE IF NOT EXISTS odb10v1_level2species (
        top_taxid integer NOT NULL,
        odb_org_id integer NOT NULL,
        hops integer NOT NULL,
        intermediate_levels text NOT NULL);
    """
    create_table(conn, sql)

    # Add the data.
    sql = """
      INSERT INTO odb10v1_level2species(top_taxid, odb_org_id, hops,
        intermediate_levels)
      VALUES(?,?,?,?);
    """
    fh = open(orthodb_dir + "/odb10v1_level2species.tab")
    for line in fh:
        cols = line.strip("\n").split("\t")
        cur = conn.cursor()
        cur.execute(sql, cols)

    # Create the indexes.
    sql = "CREATE INDEX odb10v1_level2species_idx1 ON odb10v1_level2species(odb_org_id);"
    cur.execute(sql)

    # Commit the transaction.
    conn.commit()


def index_odb10v1_species(conn, orthodb_dir):
    """Creates and inserts data into the odb10v1_species table.
    :param conn: Connection object
    :param orthodb_dir: The directory where the orthodb files are kept.
    :return:
    """
    print("Indexing odb10v1_species.tab")

    # Create the species table.
    sql = """
      CREATE TABLE IF NOT EXISTS odb10v1_species (
        ncbi_taxid integer NOT NULL,
        odb_org_id integer NOT NULL,
        scientific_name text NOT NULL,
        genome_assembly_id text,
        num_clustered_genes integer NOT NULL,
        num_OGs integer NOT NULL,
        mapping_type text NOT NULL);
    """
    create_table(conn, sql)

    # Add the data.
    sql = """
      INSERT INTO odb10v1_species(ncbi_taxid, odb_org_id, scientific_name,
        genome_assembly_id, num_clustered_genes, num_OGs, mapping_type)
      VALUES(?,?,?,?,?,?,?);
    """
    fh = open(orthodb_dir + "/odb10v1_species.tab")
    for line in fh:
        cols = line.strip("\n").split("\t")
        cur = conn.cursor()
        cur.execute(sql, cols)
    conn.commit()

    # Create the indexes.
    sql = "CREATE INDEX odb10v1_species_idx1 ON odb10v1_species(ncbi_taxid);"
    cur.execute(sql)
    sql = "CREATE INDEX odb10v1_species_idx2 ON odb10v1_species(odb_org_id);"
    cur.execute(sql)

    # Commit the transaction.
    conn.commit()


if __name__ == "__main__":
    """
    The main subroutine of the program.  The script receives a single arugment:
    the directory where the orthodb files are kept.
    """

    script_name, orthodb_dir = sys.argv

    if not os.path.exists(orthodb_dir):
        print("ERROR: The orthodb data directory is not present at {}".format(orthodb_dir))
        exit(1)
    if not os.path.exists(orthodb_dir + "/odb10v1_species.tab"):
        print("ERROR: The odb10v1_species.tab file is not present")
        exit(1)
    if not os.path.exists(orthodb_dir + "/odb10v1_level2species.tab"):
        print("ERROR: The odb10v1_level2species.tab file is not present")
        exit(1)
    if not os.path.exists(orthodb_dir + "/odb10v1_OGs.tab"):
        print("ERROR: The odb10v1_OGs.tab file is not present")
        exit(1)
    if not os.path.exists(orthodb_dir + "/odb10v1_OG2genes.tab"):
        print("ERROR: The odb10v1_OG2genes.tab file is not present")
        exit(1)
    if not os.path.exists(orthodb_dir + "/odb10v1_OG_xrefs.tab"):
        print("ERROR: The odb10v1_OG_xrefs.tab file is not present")
        exit(1)

    # Cleanup any existing database. We'll start over if this script is called
    # after it has been run repviously.
    if os.path.exists(orthodb_dir + "/odb10v1.db"):
        os.remove(orthodb_dir + "/odb10v1.db")

    # Connect to the SQLite database.
    conn = create_connection(orthodb_dir + "/odb10v1.db")

    # Index each of the data files.
    index_odb10v1_species(conn, orthodb_dir)
    index_odb10v1_level2species(conn, orthodb_dir)
    index_odb10v1_OGs(conn, orthodb_dir)
    index_odb10v1_OG2genes(conn, orthodb_dir)
    index_odb10v1_OG_xrefs(conn, orthodb_dir)

    # Close the connection.
    conn.close()

    # Return sucessful finish.
    exit(0)
