#!/usr/bin/env bash

# Prepare the directory
mkdir -p uniprot_trembl
rm -rf uniprot_tremble/*

# Download the file
cd uniprot_trembl
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_trembl.fasta.gz
gunzip uniprot_trembl.fasta.gz

# Index the file for BLAST
docker run -v ${PWD}:/Annotater/data -u $(id -u ${USER}):$(id -g ${USER}) annotater/diamond:0.9.25-0.9 /bin/bash -c "cd /Annotater/data; diamond makedb --threads 4 --in uniprot_trembl.fasta -d uniprot_trembl"
