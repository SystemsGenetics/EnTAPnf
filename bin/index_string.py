import argparse
import sqlite3


CSV_COL_SIZE = 16
ROWS_PER_DB = 620000000
DB_EXT = ".db"




def makeDB(path):
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




def populateDB(conn,ifile):
    count = 0
    while True:
        line = ifile.readline()
        if not line:
            conn.commit()
            return False
        parts = line.split()
        if len(parts)==CSV_COL_SIZE:
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
                conn.cursor().execute(
                    "CREATE INDEX idx_protein_links_full_protein1"
                    " ON protein_links_full (protein1);"
                )
                conn.cursor().execute(
                    "CREATE INDEX idx_protein_links_full_protein2"
                    " ON protein_links_full (protein2);"
                )
                conn.commit()
                return True




parser = argparse.ArgumentParser()
parser.add_argument('--links',dest='linksPath',type=str,required=True)
parser.add_argument('--out',dest='outPath',type=str,required=True)
args = parser.parse_args()
dbIdx = 0
conn = makeDB("%s%02d%s"%(args.outPath,dbIdx,DB_EXT))
with open(args.linksPath,"r") as ifile:
    ifile.readline()
    while populateDB(conn,ifile):
        dbIdx += 1
        conn = makeDB("%s%02d%s"%(args.outPath,dbIdx,DB_EXT))
