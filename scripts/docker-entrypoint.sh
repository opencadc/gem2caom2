#!/bin/bash

if [[ ! -e ${PWD}/config.yml ]]; then
  if [[ "${@}" == "gem_run_query" ]]; then
    echo "visit"
    cp /config_with_visit.yml ${PWD}/config.yml
  else
    echo "ingest"
    cp /config_with_ingest.yml ${PWD}/config.yml
  fi
fi

if [[ ! -e ${PWD}/state.yml ]]; then
  cp /state.yml ${PWD}
fi

exec "${@}"
