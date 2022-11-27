#!/usr/bin/env bash

version=`cat ../VERSION  | perl -p -e 's/^version\s*=\s*(.*)$/\1/'`

# Prepare the directory
mkdir -p string
rm -rf ./string/*

# Download the files
cd string
wget https://stringdb-static.org/download/protein.sequences.v11.0.fa.gz
gunzip protein.sequences.v11.0.fa.gz

wget https://stringdb-static.org/download/protein.links.full.v11.0.txt.gz
gunzip protein.links.full.v11.0.txt.gz

wget https://stringdb-static.org/download/protein.info.v11.0.txt.gz
gunzip protein.info.v11.0.txt.gz

# Index the file for BLAST
docker run \
  -v ${PWD}:/EnTAP/data \
  -u $(id -u ${USER}):$(id -g ${USER}) \
  quay.io/biocontainers/diamond:2.0.15--hb97b32f_0 \
   /bin/bash -c "cd /EnTAP/data; diamond makedb --threads 4 --in protein.sequences.v11.0.fa -d protein.sequences.v11.0"
