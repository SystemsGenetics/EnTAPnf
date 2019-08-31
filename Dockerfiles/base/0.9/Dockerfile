# Base image.
FROM ubuntu:16.04


################## METADATA ######################
LABEL base_image="ubuntu:16.04"
LABEL software="AnnoTater Base Image"
LABEL software.version="1.0.0"
LABEL about.summary="Base image for AnnoTater"
LABEL about.home=""
LABEL about.documentation=""
LABEL about.license_file=""

# Set noninteractive mode for apt-get.
ENV DEBIAN_FRONTEND noninteractive


# Update and install basic packages.
RUN apt-get update && apt-get upgrade -y wget curl unzip \
  && apt-get install -y git


# Tini for signal processing and zombie killing.
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini


# Create a non-root user that will run the called programs.
RUN groupadd -g 61000 bioinfo \
  && useradd -g 61000 --no-log-init --create-home --shell /bin/bash -u 61000 bioinfo

# Create a directory for annotater specif file that other
# dependent images may want to use.
RUN mkdir /annotater  && chown bioinfo /annotater

# Change the active user.
USER bioinfo

WORKDIR /home/bioinfo

# Define the command and parameters that will be executed when this
# container is first run.
ENTRYPOINT ["/tini", "--"]
