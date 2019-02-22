#!/bin/bash
for ii in $(ls -d *)
do
    if [[ -d ${ii} && ${ii} != 'votable' ]]
    then
	echo "    # $ii"
#	mkdir -p ${ii}/program
#	for jj in $(grep program_id ${ii}/json/* | awk -F'"' '{print $4}')
#	do
#		echo $jj
#		curl "https://archive.gemini.edu/programinfo/${jj}" > ${ii}/program/${jj}.xml
#	done
	for jj in $(ls ${ii}/*.in.xml)
	do
		obs_id=$( grep observationID ${jj} | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' )
		file_id=$( grep uri ${jj} | grep fits | awk -F'>' '{print $2}' | awk -F'<' '{print $1}' | awk -F'/' '{print $2}' | awk -F'.' '{print $1"."$2}' )
		program_id=$(grep program_id ${ii}/json/${obs_id}.json | awk -F'"' '{print $4}')
		echo "    '${file_id}': ['${obs_id}', '${ii}', '${program_id}'], "
	done
    fi
done
date
exit 0
