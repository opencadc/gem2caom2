#!/bin/bash

# CONTAINER="bucket.canfar.net/gem2caom2"
CONTAINER="gem_run_test"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get the container"
docker pull ${CONTAINER} || exit $?

echo "Run gem_run container"
docker run --rm --name gem_run_query -v ${PWD}:/usr/src/app/ ${CONTAINER} gem_run_query || exit $?

date
exit 0
