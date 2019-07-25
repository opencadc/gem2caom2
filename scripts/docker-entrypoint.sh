#!/bin/bash

if [[ ! -e ${PWD}/config.yml ]]; then
  if [[ "${@}" == "gem_run_state" ]]; then
    echo "ingest"
    cp /config_with_ingest.yml ${PWD}/config.yml
  else
    # the config required for gem_run_public, gem_run_query
    echo "visit"
    cp /config_with_visit.yml ${PWD}/config.yml
  fi
fi

if [[ ! -e ${PWD}/state.yml ]]; then
  if [[ "${@}" == "gem_run_public" ]]; then
    yesterday=$(date -d yesterday "+%d-%b-%Y %H:%M")
    echo "bookmarks:
    vlass_timestamp:
      last_record: $yesterday
" > ${PWD}/state.yml
  else
    cp /state.yml ${PWD}
  fi
fi

exec "${@}"
