FROM annotater/base:0.9

USER root

# Install Oracle Java 8 which is used by InterProScan. Rather than use
# the apt-add-repository program (which installs a lot of packages) we will
# manually add the repository.
RUN echo "deb http://ppa.launchpad.net/openjdk-r/ppa/ubuntu xenial main" >> /etc/apt/sources.list \
  && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys 86F44E2A

RUN apt-get update \
  && apt-get -y install openjdk-8-jdk \
  && update-alternatives --config java

# Install Pyhon3
RUN apt-get install -qq -y python3-pip \
  && pip3 install -q numpy pandas xmltodict

# Install InterProScan
RUN wget ftp://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.36-75.0/alt/interproscan-core-5.36-75.0.tar.gz \
  && tar -zxvf interproscan-core-5.36-75.0.tar.gz \
  && mv interproscan-5.36-75.0 /usr/local/interproscan

# Install IPRScan dependency
RUN apt-get install -y elfutils

USER bioinfo
