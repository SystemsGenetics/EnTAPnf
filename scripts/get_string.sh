#!/usr/bin/env bash

# Prepare the directory
mkdir -p string
rm -rf ./string/*

# Download the files
cd string
wget https://stringdb-static.org/download/protein.sequences.v11.0.fa.gz


# Index the file for BLAST
docker run -v ${PWD}:/Annotater/data -u $(id -u ${USER}):$(id -g ${USER}) annotater/diamond:0.9.25-0.9 /bin/bash -c "cd /Annotater/data; diamond makedb --threads 4 --in protein.sequences.v11.0.fa -d protein.sequences.v11.0"
