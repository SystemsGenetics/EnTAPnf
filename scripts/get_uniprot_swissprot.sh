#!/usr/bin/env bash

# Prepare the directory
mkdir -p uniprot_sprot
rm -rf uniprot_sprot/*
cd uniprot_sprot

# Download the FASTA file
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Download the protein ID to EC number mapping file
wget ftp://ftp.expasy.org/databases/enzyme/enzyme.dat

# Index the file
docker run \
  -v ${PWD}:/EnTAP/data \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  systemsgenetics/entap:flask \
   /bin/bash -c "cd /EnTAP/data/; diamond makedb --threads 4 --in uniprot_sprot.fasta -d uniprot_sprot"
