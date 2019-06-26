#!/bin/bash

CONTAINER="bucket.canfar.net/gem2caom2"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Get the container"
docker pull ${CONTAINER} || exit $?

echo "Run gem_run container"
docker run -m=7g --rm --name gem_run_query -v ${PWD}:/usr/src/app/ ${CONTAINER} gem_run_state || exit $?

date
exit 0
