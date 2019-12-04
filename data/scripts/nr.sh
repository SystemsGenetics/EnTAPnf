#!/usr/bin/env bash

# Prepare the directory
mkdir -p nr
rm -rf nr/*
cd nr

# Download the file
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr

# Index the file for BLAST
docker run -v ${PWD}:/Annotater/data -u $(id -u ${USER}):$(id -g ${USER}) annotater/diamond:0.9.25-0.9 /bin/bash -c "cd /Annotater/data; diamond makedb --threads 4 --in nr -d nr"
