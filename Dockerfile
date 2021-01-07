FROM opencadc/pandas:3.8-slim

RUN apt-get update -y && apt-get dist-upgrade -y && \
    apt-get install -y build-essential \
                       git && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN pip install bs4 \
    cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    caom2utils \
    deprecated \
    ftputil \
    importlib-metadata \
    matplotlib \
    pillow \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

RUN mkdir /app && mkdir /app/data

ADD https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/cadcsw/2019-07-03_from_paul.txt /app/data/from_paul.txt

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git --branch ${OPENCADC_BRANCH} --single-branch && \
    rm -rf ./caom2tools/.git && \
    rm -rf ./caom2tools/caom2 && \
    rm -rf ./caom2tools/caom2repo && \
    rm -rf ./caom2tools/caom2utils/caom2utils/tests && \
    pip install ./caom2tools/caom2utils

RUN git clone https://github.com/${OPENCADC_REPO}/caom2pipe.git --branch ${OPENCADC_BRANCH} --single-branch && \
    rm -rf ./caom2pipe/.git && \
    rm -rf ./caom2pipe/caom2pipe/tests && \
    pip install ./caom2pipe

RUN git clone https://github.com/${OPENCADC_REPO}/gem2caom2.git --branch ${OPENCADC_BRANCH} --single-branch && \
    rm -rf ./gem2caom2/.git && \
    rm -rf ./gem2caom2/gem2caom2/test && \
    pip install ./gem2caom2 && \
    cp ./gem2caom2/scripts/docker-entrypoint.sh / && \
    cp ./gem2caom2/scripts/config_with_ingest.yml / && \
    cp ./gem2caom2/scripts/config_with_visit.yml / && \
    cp ./gem2caom2/scripts/state.yml /

ENTRYPOINT ["/docker-entrypoint.sh"]
