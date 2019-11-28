#!/bin/bash

COLLECTION="gem"
IMAGE="bucket.canfar.net/${COLLECTION}2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run -m=7g --rm --name ${COLLECTION}_validator -v ${PWD}:/usr/src/app/ ${IMAGE} ${COLLECTION}_validate || exit $?

date
exit 0
