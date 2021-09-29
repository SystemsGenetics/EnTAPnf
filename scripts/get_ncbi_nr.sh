#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir -p nr
rm -rf nr/*
cd nr

# Download the file
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr

# Index the file for BLAST
docker run -v ${PWD}:/EnTAP/data -u $(id -u ${USER}):$(id -g ${USER}) annotater/diamond:0.9.25-${version} /bin/bash -c "cd /EnTAP/data; diamond makedb --threads 4 --in nr -d nr"
