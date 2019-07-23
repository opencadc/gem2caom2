#!/bin/bash

CONTAINER="bucket.canfar.net/gem2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get the container"
docker pull ${CONTAINER} || exit $?

echo "Run container ${container}"
docker run -m=7g --rm --name gem_run_public -v ${PWD}:/usr/src/app/ ${CONTAINER} gem_run_public || exit $?

date
exit 0
