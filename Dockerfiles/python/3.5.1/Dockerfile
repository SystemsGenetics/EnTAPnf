FROM annotater/base:0.9

USER root

RUN apt-get install -qq -y python3-pip \
  && pip3 install -q numpy pandas xmltodict

RUN pip3 install -q biopython

USER bioinfo
