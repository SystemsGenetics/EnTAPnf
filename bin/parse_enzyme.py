import sys
from parse_blastxml import parseBLASTXMLfile
from _enzymelookup import EnzymeLookup




data = parseBLASTXMLfile(sys.argv[1],True)
query_def = data["Query_def"]
hit_id = data["Hit_id"]
enzymeLookup = EnzymeLookup(sys.argv[2])
for (query,hit) in zip(query_def,hit_id):
    hit_parts = hit.split("|")
    enzymeLookup.find(query.split()[0],hit_parts[1],hit_parts[2])
