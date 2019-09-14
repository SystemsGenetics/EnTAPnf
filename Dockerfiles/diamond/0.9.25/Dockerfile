FROM systemsgenetics/annotater:0.9

USER root

# Install NCBI Blast+
RUN wget https://github.com/bbuchfink/diamond/releases/download/v0.9.25/diamond-linux64.tar.gz \
  && tar -zxvf diamond-linux64.tar.gz \
  && mv diamond /usr/local/bin/diamond


USER bioinfo
