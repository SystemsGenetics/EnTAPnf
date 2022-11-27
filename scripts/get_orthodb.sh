#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir orthodb
rm -rf ./orthodb/*
cd orthodb

# Get the data
wget https://v101.orthodb.org/download/odb10v1_level2species.tab.gz
gunzip odb10v0_level2species.tab.gz

wget https://v101.orthodb.org/download/odb10v1_OG2genes.tab.gz
gunzip odb10v0_OG2genes.tab.gz

wget https://v101.orthodb.org/download/odb10v1_OGs.tab.gz
gunzip odb10v0_OGs.tab.gz

wget https://v101.orthodb.org/download/odb10v1_all_og_fasta.tab.gz
gunzip odb10_all_og_fasta.tab.gz

wget https://v101.orthodb.org/download/odb10v1_OG_xrefs.tab.gz
gunzip odb10v1_OG_xrefs.tab.gz

wget https://v101.orthodb.org/download/odb10v1_species.tab.gz
gunzip odb10v1_species.tab.gz

wget https://v101.orthodb.org/download/odb10v1_gene_xrefs.tab.gz
gunzip odb10v1_gene_xrefs.tab.gz

# Index the file for BLAST
docker run \
  -v ${PWD}:/EnTAP/data \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  quay.io/biocontainers/diamond:2.0.15--hb97b32f_0 \
  /bin/bash -c "cd /EnTAP/data; diamond makedb --threads 4 --in odb10v1_all_og_fasta.tab -d odb10v1_all_og"

docker run \
  -v ${PWD}:/EnTAP/data \
  -v ${PWD}/../../bin:/EnTAP/bin \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  annotater/python:3.7-${version} \
  /bin/bash -c "cd /EnTAP/data; /EnTAP/bin/index_orthodb.py ."
