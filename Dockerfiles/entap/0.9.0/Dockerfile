FROM annotater/base:0.9

USER root

# Install gcc and cmake dependencies
RUN apt-get -y install build-essential && apt-get -y install cmake

# Install Python dependencies
RUN apt-get -y install python3 python3-pip && python3 -m pip install -U pip \
  && python3 -m pip install -U matplotlib


USER bioinfo

RUN git clone https://gitlab.com/enTAP/EnTAP.git && cd EnTap \
  && git checkout v0.9.0-beta
