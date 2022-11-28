#!/bin/bash

outfile=$1

# Get the list of blast txt files
txt_files=`ls *.txt`

# Write output compatible with this blast format option:
# qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovhsp stitle
touch $outfile
for f in $txt_files; do
    cat $f | awk -F"\t" 'BEGIN {OFS="\t"} {print $1,$4,$12,$19,$18,$15,$20,$21,$22,$23,$10,$8,$17,$5}'| grep -v Query_id | grep -v -P "^[\s]+\$" >> $outfile
done;

exit 0
