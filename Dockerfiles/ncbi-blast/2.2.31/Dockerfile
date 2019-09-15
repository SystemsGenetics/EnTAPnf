FROM annotater/base:0.9

USER root

# Install NCBI Blast+
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.2.31/ncbi-blast-2.2.31+-x64-linux.tar.gz \
  && tar -zxvf ncbi-blast-2.2.31+-x64-linux.tar.gz \
  && mv ncbi-blast-2.2.31+ /usr/local/ncbi-blast


USER bioinfo
