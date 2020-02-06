#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir orthodb
rm -rf ./orthodb/*
cd orthodb

# Get the data
wget https://v100.orthodb.org/download/odb10v0_level2species.tab.gz
gunzip odb10v0_level2species.tab.gz

wget https://v100.orthodb.org/download/odb10v0_OG2genes.tab.gz
gunzip odb10v0_OG2genes.tab.gz

wget https://v100.orthodb.org/download/odb10v0_OGs.tab.gz
gunzip odb10v0_OGs.tab.gz

wget https://v100.orthodb.org/download/odb10_all_og_fasta.tab.gz
gunzip odb10_all_og_fasta.tab.gz

# Index the file for BLAST
docker run -v ${PWD}:/Annotater/data -u $(id -u ${USER}):$(id -g ${USER}) annotater/diamond:0.9.25-${version} /bin/bash -c "cd /Annotater/data; diamond makedb --threads 4 --in odb10_all_og_fasta.tab -d odb10_all_og"
