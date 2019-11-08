# Prepare the directory
mkdir -p nr
rm -rf nr/*

# Download the file
cd nr
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr

# Index the file for BLAST
diamond makedb --threads 4 --in nr.faa -d nr
