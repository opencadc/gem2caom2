#!/bin/bash

curl https://archive.gemini.edu/jsonsummary/canonical/filepre=${1} > ${1}.json
exit 0
