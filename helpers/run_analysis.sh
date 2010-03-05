#!/bin/bash
#
HGSC_BF_PROD_PIPE="/stornext/snfs1/next-gen/solid/hgsc.solid.pipeline"

#rm -rf lsf_logs/ output/ cluster_JOBS.sh reads metric* track* rg.txt marked*.txt
mkdir ./reads
$HGSC_BF_PROD_PIPE/bin/bfast.split.reads.rb ./bf.config.yaml
$HGSC_BF_PROD_PIPE/bin/bfast.lsf.submit.rb ./bf.config.yaml
exit 

bkill 0
echo "Sending pipeline jobs to cluster ..."
./cluster_JOBS.sh &> /dev/null
