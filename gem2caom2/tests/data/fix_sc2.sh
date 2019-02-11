#!/bin/bash

for ii in GMOS/GN-CAL20150925-2-007.actual.xml GMOS/GN-CAL20150217-2-003.actual.xml GMOS/GN-2015A-Q-91-5-002.actual.xml GMOS/GN-2015A-C-4-24-086.actual.xml GMOS/GN-2015A-C-2-96-002.actual.xml GMOS/GN-2013B-Q-28-150-002.actual.xml
do
        obs_id=$(echo $ii | awk -F'/' '{print $2}' | awk -F '.' '{print $1}')
        caom2-repo delete --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo GEMINI ${obs_id}
        caom2-repo create --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo ${ii}
done
