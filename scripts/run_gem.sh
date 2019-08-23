#!/bin/bash
if [[ $# -ne 1 ]]
then
	echo "arg is obs id"
	exit 1
fi

obs_id=${1}
fid=$(grep ${obs_id} /app/data/from_paul.txt | awk -F\| '{print $2}' | awk -F'.fits' '{print $1}' | xargs)
echo "${obs_id} ${fid}"

module="/usr/src/app/gem2caom2/gem2caom2/main_app.py"
external="--external_url "
lineage="--lineage "
for ii in ${fid}; do
  external="${external} https://archive.gemini.edu/fullheader/${ii}.fits"
  lineage="${lineage} ${ii}/gemini:GEM/${ii}.fits"
done

# echo "gem2caom2 --no_validate --observation GEMINI ${obs_id} --external_url https://archive.gemini.edu/fullheader/${fid}.fits --plugin /app/gem2caom2/gem2caom2/main_app.py --module /app/gem2caom2/gem2caom2/main_app.py --out ./${fid}.xml --lineage ${fid}/gemini:GEM/${fid}.fits"
gem2caom2 --no_validate --observation GEMINI ${obs_id} ${external} --plugin ${module} --module ${module} --out ./${obs_id}.xml ${lineage}
