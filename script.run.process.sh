#!/bin/bash
#
# Run the bfast jobs via this script. It will create ok/error files
# that we can use later to figure out what jobs complited and which
# didn't
#

help()
{
	echo "Usage:"
	echo "$0 job_name cmd"
	exit 1
}

[ ".$1" == "." ] && help
[ ".$2" == "." ] && help

job_name="$1"
cmd="$2"
output_dir="./track_jobs"

mkdir -p $output_dir
$cmd
if [ $? -eq 0 ]
then
	touch $output_dir/$job_name.ok
else
	touch $output_dir/$job_name.error
fi
