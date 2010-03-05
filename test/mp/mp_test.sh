#!/bin/bash
# 
set -e

rm -rf lsf_logs/ output/ cluster_JOBS.sh reads metric* track* rg.txt marked*.txt
mkdir ./reads
../../bin/bfast.split.reads.rb ./bf.config.yaml
../../bin/bfast.lsf.submit.rb ./bf.config.yaml
exit # !!!!!!!!!!!


bkill 0
echo "Sending pipeline jobs to cluster ..."
./cluster_JOBS.sh &> /dev/null
