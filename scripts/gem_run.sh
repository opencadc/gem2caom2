#!/bin/bash

IMAGE="bucket.canfar.net/gem2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE}

echo "Run image ${IMAGE}"
docker run --rm --name gem_run_builder -v ${PWD}:/usr/src/app/ ${IMAGE} gem_run || exit $?

date
exit 0
