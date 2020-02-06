#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir -p uniprot_sprot
rm -rf uniprot_sprot/*
cd uniprot_sprot

# Download the file
wget ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
gunzip uniprot_sprot.fasta.gz

# Index the file
docker run -v ${PWD}:/Annotater/data -u $(id -u ${USER}):$(id -g ${USER}) annotater/diamond:0.9.25-${version} /bin/bash -c "cd /Annotater/data; diamond makedb --threads 4 --in uniprot_sprot.fasta -d uniprot_sprot"
