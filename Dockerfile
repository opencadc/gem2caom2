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
    importlib-metadata \
    matplotlib \
    pillow \
    python-dateutil \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

RUN mkdir /app && mkdir /app/data

ADD https://www.cadc-ccda.hia-iha.nrc-cnrc.gc.ca/files/vault/cadcsw/2019-07-03_from_paul.txt /app/data/from_paul.txt

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/gem2caom2@${PIPE_BRANCH}#egg=gem2caom2

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
