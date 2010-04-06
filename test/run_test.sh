#!/bin/bash

# This script is the master script to test three different analysis types
# namely fragment (fr), mate-pair (mp) and capture (cap).

test_type=""

help()
{
  echo "usage:"
  echo "$0 <TestType>"
  echo "    <TestType> - Type of run. Allowed values :"
  echo "                 fr, mp, cap"
  exit 1
}

[ $# -ne 1 ] && help
test_type=`echo $1 | tr [:upper:] [:lower:]`

[ "$test_type" != "mp" ] && [ "$test_type" != "fr" ] &&
 [ "$test_type" != "cap" ] && help

rm -rf $test_type

mkdir -p $test_type
cp generate_yaml.rb $test_type

cd $test_type

mkdir ./reads
mkdir ./track_jobs
echo "Generating YAML"

ruby19 "generate_yaml.rb" $test_type 

../../bin/bfast.split.reads.rb ./bf.config.yaml
../../bin/bfast.lsf.submit.rb ./bf.config.yaml

echo "Sending pipeline jobs to cluster ..."
./cluster_JOBS.sh
cd ../
