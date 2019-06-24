FROM python:3.6-alpine

RUN apk --no-cache add \
        bash \
        coreutils \
        gcc \
        git \
        g++ \
        libffi-dev \
        libmagic \
        libxml2-dev \
        libxslt-dev \
        make \
        musl-dev \
        openssl-dev

RUN pip install aenum && \
    pip install astropy && \
    pip install cadcdata && \
    pip install cadctap && \
    pip install caom2repo && \
    pip install funcsigs && \
    pip install future && \
    pip install numpy && \
    pip install PyYAML && \
    pip install spherical-geometry && \
    pip install xml-compare

WORKDIR /usr/src/app

RUN pip install bs4

RUN apk --no-cache add imagemagick

RUN git clone https://github.com/SharonGoliath/caom2tools.git && \
  cd caom2tools && git pull origin master && \
  pip install ./caom2 && \
  pip install ./caom2utils && pip install ./caom2pipe && cd ..

RUN git clone https://github.com/SharonGoliath/gem2caom2.git && \
  pip install ./gem2caom2 && \
  cp ./gem2caom2/scripts/docker-entrypoint.sh / && \
  cp ./gem2caom2/scripts/config_with_ingest.yml / && \
  cp ./gem2caom2/scripts/config_with_visit.yml / && \
  cp ./gem2caom2/scripts/state.yml /

RUN mkdir /app && mkdir /app/data

COPY ./2018-12-17_from_paul.txt /app/data/from_paul.txt

RUN apk --no-cache del git

ENTRYPOINT ["/docker-entrypoint.sh"]

