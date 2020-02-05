#!/usr/bin/env bash

# Preapre the data directory.
mkdir -p inteproscan
rm -rf interproscan/*
cd interproscan

# Get the InterProScan application nad data.
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/interproscan-5.36-75.0-64-bit.tar.gz
tar -zxvf interproscan-5.36-75.0-64-bit.tar.gz

# Add the Panter database.
cd data
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-14.1.tar.gz
tar -zxvf panther-data-14.1.tar.gz
