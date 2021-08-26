#!/bin/bash

# IMAGE="bucket.canfar.net/gem2caom2"
IMAGE="gem_build_test"

echo "Get a proxy certificate"
cp $HOME/.ssl/cadcproxy.pem ./ || exit $?

echo "Pull ${IMAGE}"
# docker pull ${IMAGE} || exit $?

echo "Run ${IMAGE}"
docker run --rm --name gem_run_edu_query -v ${PWD}:/usr/src/app/ ${IMAGE} gem_run_edu_query || exit $?

date
exit 0
