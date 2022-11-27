#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir -p uniprot_trembl
rm -rf uniprot_tremble/*

# Download the file
cd uniprot_trembl
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz

# Index the file for BLAST
docker run \
  -v ${PWD}:/EnTAP/data \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  quay.io/biocontainers/diamond:2.0.15--hb97b32f_0 \
  /bin/bash -c "cd /EnTAP/data; diamond makedb --threads 4 --in uniprot_trembl.fasta -d uniprot_trembl"
