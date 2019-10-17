# Base image.
FROM ubuntu:18.04


################## METADATA ######################
LABEL base_image="ubuntu:18.04"
LABEL software="Annotater Base Image"
LABEL software.version="1.0.0"
LABEL about.summary="Base Docker image for Annotater"
LABEL about.home=""
LABEL about.documentation=""
LABEL about.license_file=""


# Set noninteractive mode for apt-get.
ENV DEBIAN_FRONTEND noninteractive


# Update and install basic packages.
RUN apt-get update -qq \
  && apt-get install -qq -y curl git unzip wget gnupg


# Tini for signal processing and zombie killing.
ENV TINI_VERSION v0.18.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini


# Define the command and parameters that will be executed when this
# container is first run.
ENTRYPOINT ["/tini", "--"]
