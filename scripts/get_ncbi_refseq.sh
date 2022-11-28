#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir -p refseq/plant
rm -rf refseq/plant*
cd refseq

# Download the amino acid sequence files
wget -r -A '*.protein.faa.gz' ftp://ftp.ncbi.nlm.nih.gov/refseq/release/plant/
gunzip ./ftp.ncbi.nlm.nih.gov/refseq/release/plant/*.gz

# Combine the protein sequences into one large file.
cat ./ftp.ncbi.nlm.nih.gov/refseq/release/plant/*.faa > refseq_plant.protein.faa

# Remove the download directory.
rm -rf ./ftp.ncbi.nlm.nih.gov/

# Index the file for Diamond BLAST
docker run \
  -v ${PWD}:/EnTAP/data \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  systemsgenetics/entap:flask \
  /bin/bash -c "cd /EnTAP/data; diamond makedb --threads 4 --in refseq_plant.protein.faa -d refseq_plant.protein"
