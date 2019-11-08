mkdir -p inteproscan
rm -rf interproscan/*
cd interproscan
wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/interproscan-5.36-75.0-64-bit.tar.gz
tar -zxvf interproscan-5.36-75.0-64-bit.tar.gz

#wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/data/panther-data-12.0.tar.gz
tar -zxvf panther-data-12.0.tar.gz

