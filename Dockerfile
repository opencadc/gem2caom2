FROM opencadc/pandas:3.10-slim as builder

RUN apt-get update --no-install-recommends  && apt-get dist-upgrade -y && \
    apt-get install -y build-essential \
                       libcfitsio-dev \
                       git \
                       imagemagick \
                       libcurl4-openssl-dev \
                       libnsl-dev \
                       libssl-dev \
                       zlib1g-dev && \
    rm -rf /var/lib/apt/lists/ /tmp/* /var/tmp/*

ARG OPENCADC_BRANCH=master
ARG OPENCADC_REPO=opencadc
ARG FITSVERIFY_VERSION=4.20
ARG FITSVERIFY_URL=https://heasarc.gsfc.nasa.gov/docs/software/ftools/fitsverify/fitsverify-${FITSVERIFY_VERSION}.tar.gz

ADD ${FITSVERIFY_URL} /usr/local/src

RUN cd /usr/local/src \
  && tar axvf fitsverify-${FITSVERIFY_VERSION}.tar.gz \
  && cd fitsverify-${FITSVERIFY_VERSION} \
  && gcc -o fitsverify ftverify.c fvrf_data.c fvrf_file.c fvrf_head.c fvrf_key.c fvrf_misc.c -DSTANDALONE -I/usr/local/include -L/usr/local/lib -lcfitsio -lm -lnsl \
  && cp ./fitsverify /usr/local/bin/ \
  && ldconfig

WORKDIR /usr/src/app

RUN git clone https://github.com/opencadc/cadctools.git && \
    cd cadctools && \
    pip install ./cadcutils && \
    pip install ./cadcdata && \
    cd ..

RUN git clone https://github.com/${OPENCADC_REPO}/caom2tools.git && \
    cd caom2tools && \
    git checkout ${OPENCADC_BRANCH} && \
    pip install ./caom2utils && \
    cd ..

RUN pip install git+https://github.com/${OPENCADC_REPO}/caom2pipe@${OPENCADC_BRANCH}#egg=caom2pipe

RUN pip install git+https://github.com/${OPENCADC_REPO}/gem2caom2@${OPENCADC_BRANCH}#egg=gem2caom2

FROM python:3.10-slim
WORKDIR /usr/src/app

COPY --from=builder /usr/local/lib/python3.10/site-packages/ /usr/local/lib/python3.10/site-packages/
COPY --from=builder /usr/local/bin/* /usr/local/bin/
COPY --from=builder /usr/local/.config/* /usr/local/.config/

COPY --from=builder /etc/magic /etc/magic
COPY --from=builder /etc/magic.mime /etc/magic.mime
COPY --from=builder /usr/lib/x86_64-linux-gnu/libmagic* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/file/magic.mgc /usr/lib/file/
COPY --from=builder /usr/share/misc/magic /usr/share/misc/magic
COPY --from=builder /usr/share/misc/magic.mgc /usr/share/misc/magic.mgc
COPY --from=builder /usr/share/file/magic.mgc /usr/share/file/magic.mgc

# fitsverify
COPY --from=builder /usr/lib/x86_64-linux-gnu/libcfitsio* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libcurl-gnutls* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libnghttp2* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/librtmp* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libssh2* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libpsl* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libldap* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/liblber* /usr/lib/x86_64-linux-gnu/
COPY --from=builder /usr/lib/x86_64-linux-gnu/libsasl* /usr/lib/x86_64-linux-gnu/

RUN useradd --create-home --shell /bin/bash cadcops
RUN chown -R cadcops:cadcops /usr/src/app
USER cadcops

ENTRYPOINT ["/usr/local/bin/docker-entrypoint.sh"]
