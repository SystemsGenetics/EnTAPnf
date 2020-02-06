#!/usr/bin/env bash

./nextflow main.nf -profile standard -with-docker -resume

# To force a full re-run of the  test comment the line above and uncomment the following
#rm -rf work .nextflow .nextflow.log* ./output
#./nextflow main.nf -profile standard -with-docker
