FROM opencadc/pandas:3.8-slim

RUN pip install cadcdata && \
    pip install cadctap && \
    pip install caom2 && \
    pip install caom2repo && \
    pip install caom2utils && \
    pip install deprecated && \
    pip install ftputil && \
    pip install PyYAML && \
    pip install pytz && \
    pip install spherical-geometry && \
    pip install vos

WORKDIR /usr/src/app

RUN apt-get update -y && apt-get dist-upgrade -y

RUN apt-get install -y git

RUN pip install bs4

RUN apt-get install -y imagemagick

RUN mkdir /app && mkdir /app/data

ADD https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/cadcsw/2019-07-03_from_paul.txt /app/data/from_paul.txt

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG OMC_REPO=opencadc-metadata-curation

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git --branch ${OPENCADC_BRANCH} --single-branch && \
    pip3 install ./caom2tools/caom2 && \
    pip3 install ./caom2tools/caom2utils

RUN git clone https://github.com/${OMC_REPO}/caom2pipe.git && \
  pip install ./caom2pipe

RUN git clone https://github.com/${OMC_REPO}/gem2caom2.git && \
  pip install ./gem2caom2 && \
  cp ./gem2caom2/scripts/docker-entrypoint.sh / && \
  cp ./gem2caom2/scripts/config_with_ingest.yml / && \
  cp ./gem2caom2/scripts/config_with_visit.yml / && \
  cp ./gem2caom2/scripts/state.yml /

ENTRYPOINT ["/docker-entrypoint.sh"]

