# Prepare the directory
mkdir -p uniprot_sprot
rm -rf uniprot_sprot/*

# Download the file
cd uniprot_sprot
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Index the file
diamond makedb --threads 4 --in uniprot_sprot.fasta -d uniprot_sprot
