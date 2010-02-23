#!/bin/bash
#
# Run the bfast jobs via this script. It will create ok/error files
# that we can use later to figure out what jobs complited and which
# didn't
#

help()
{
	echo "$1"
	echo "Usage:"
	echo "$0 job_name cmd"
	exit 1
}

track_dir=`cat bf.config.yaml |\
grep trackdir | awk -F: '{print $2}' | sed -e 's/ //g'`

[ ".$1" == "." ] && help "Error param: job_name"
[ ".$2" == "." ] && help "Error param: cmd"
[ ".$track_dir" == "." ] && help "Error: trackdir not found in config"

job_name="$1"
cmd="$2"
mkdir -p $track_dir

#echo "jname = $job_name "
#echo "cmd   = $cmd"

eval "$cmd"
if [ $? -eq 0 ]
then
	touch $track_dir/$job_name.ok
	exit 0
else
	touch $track_dir/$job_name.error
	exit 1
fi
