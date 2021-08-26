#!/bin/bash

IMAGE="bucket.canfar.net/gem2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get image ${IMAGE}"
docker pull ${IMAGE} || exit $?

echo "Run image ${IMAGE}"
docker run --rm --name gem_run_public -v ${PWD}:/usr/src/app/ ${IMAGE} gem_run_public || exit $?

date
exit 0
