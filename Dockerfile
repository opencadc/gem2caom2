FROM opencadc/pandas:3.8-slim

RUN apt-get update --no-install-recommends && \
    apt-get install -y build-essential git && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

RUN pip install bs4 \
    cadcdata \
    cadctap \
    caom2 \
    caom2repo \
    caom2utils \
    importlib-metadata \
    matplotlib \
    pillow \
    python-dateutil \
    PyYAML \
    spherical-geometry \
    vos

WORKDIR /usr/src/app

RUN mkdir /app && mkdir /app/data

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG PIPE_BRANCH=master
ARG PIPE_REPO=opencadc

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${PIPE_REPO}/gem2caom2@${PIPE_BRANCH}#egg=gem2caom2

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
