#!/bin/bash

echo "${@}"
script_name=$(basename $0)
if [[ ! -e ${PWD}/config.yml ]]; then
  if [[ ${script_name} == "gem_run_query" ]]; then
    cp /config_with_visit.yml ${PWD}/config.yml
  else
    cp /config_with_ingest.yml ${PWD}/config.yml
  fi
fi

if [[ ! -e ${PWD}/state.yml ]]; then
  cp /state.yml ${PWD}
fi

exec "${@}"
