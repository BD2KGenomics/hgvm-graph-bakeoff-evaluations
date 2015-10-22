#!/bin/bash
# test if vg file contains a path with given name

set -e

VG_PATH=$1
NAME=$2
NUM_PATHS=`vg view $VG_PATH -j | jq .path | jq 'map(select(.name == "'"${NAME}"'"))' | jq length`
if [ $NUM_PATHS -eq 1 ]
then
	echo 1
else
	 echo 0
fi

