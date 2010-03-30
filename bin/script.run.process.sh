#!/bin/bash
#
# Run the bfast jobs via this script. It will create ok/error files
# that we can use later to figure out what jobs complited and which
# didn't
#

mail_to="pellon@bcm.edu deiros@bcm.edu niravs@bcm.edu"

help()
{
	echo "$1"
	echo "Usage:"
	echo "$0 job_name cmd"
	exit 1
}

send_email()
{
  current_dir=$1
  job_name=$2

  (
  cat <<EOF
  pwd     : $current_dir
  job_name: $job_name
EOF
  ) | mail -s "BFAST ERROR: $job_name" $mail_to
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
  send_email `pwd` $job_name
	exit 1
fi
