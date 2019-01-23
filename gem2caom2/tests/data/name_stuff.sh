#!/bin/bash
# create the python dict that is used for looking up the relationship between obs id and file name
# in test data
for ii in GMOS
do
	for jj in $(ls ${ii}/*.in.xml)
	do
		obs_id=$( grep observationID ${jj} | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' )
		file_id=$( grep uri ${jj} | grep fits | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | awk -F'/' '{print $2}' | awk -F'.' '{print $1}' )
		echo "    '${file_id}': ['${obs_id}', 'x', '${ii}'], "
	done
done
date
exit 0
