#!/bin/bash
# 
set -e

for i in fr mp capture
do
  cd $i
  rm -rf lsf_logs/ output/ cluster_JOBS.sh reads metric* track* rg.txt marked*.txt \
  run_small* cap_stats*
  cd ..
done
