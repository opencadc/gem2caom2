#!/bin/bash

for instrument in GMOS NIRI GPI F2 GSAOI NICI
do
    echo $instrument
    for ii in $(ls ${instrument}/*.actual.xml)
    do
        echo $ii
        caom2-repo update --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo "${ii}"
    done
done
