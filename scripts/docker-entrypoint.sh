#!/bin/bash

if [[ ! -e ${PWD}/config.yml ]]; then
  if [[ "${@}" == "gem_run_state" ]]; then
    echo "ingest"
    cp /usr/local/bin/config_with_ingest.yml ${PWD}/config.yml
  else
    # the config required for gem_run_public, gem_run_query
    echo "visit"
    cp /usr/local/bin/config_with_visit.yml ${PWD}/config.yml
  fi
fi

if [[ ! -e ${PWD}/state.yml ]]; then
  if [[ "${@}" == "gem_run_public" ]]; then
    yesterday=$(date -d yesterday "+%d-%b-%Y %H:%M")
    echo "bookmarks:
    gemini_timestamp:
      last_record: $yesterday
" > ${PWD}/state.yml
  else
    cp /usr/local/bin/state.yml ${PWD}
  fi
fi

if [[ ! -e ${PWD}/cache.yml ]]; then
  cp /usr/local/bin/cache.yml ${PWD}
fi

exec "${@}"
