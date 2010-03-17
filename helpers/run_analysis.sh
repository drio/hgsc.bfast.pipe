#!/bin/bash
#
prod_pipe="/stornext/snfs1/next-gen/solid/hgsc.solid.pipeline/hgsc.bfast.pipe"
run_process="$prod_pipe/bin/script.run.process.sh"
config_file="./bf.config.yaml"
rnds=`cat /dev/urandom | tr -dc 'a-zA-Z0-9' | fold -w 9 | head -1 | cut -c1-5`
queue="normal"
lsf_logs_dir="./lsf_logs"

stuff_to_clean="lsf_logs/ output/ cluster_JOBS.sh "
stuff_to_clean="$stuff_to_clean reads metric* track* rg.txt marked*.txt"

error()
{
  echo "ERROR: $1"
  exit 1
}

to_lsf()
{
  job_name=$1
  cmd=$2
  dep_jn=$3

  if [ ".$dep_jn" == "." ]
  then
    w_dep=""
  else
    w_dep="-w 'done($dep_jn)'"
  fi

  echo "#$job_name"
  cat <<-EOF
  bsub -o $lsf_logs_dir/${job_name}.out \
  -e $lsf_logs_dir/${job_name}.err $w_dep \
  -J $job_name  \
  -q $queue  \
  "$run_process $job_name '$cmd'"
EOF
  echo ""
}

[ -f "./reads" ] && error "Reads directory exists (./reads)"
[ ! -f $config_file ] && error "Can find config file: $config"

rm -rf lsf_logs/ track* output/ cluster_JOBS.sh reads metric* track* rg.txt marked*.txt
mkdir -p ./reads
mkdir -p $lsf_logs_dir

run_name=`grep run_name $config_file | awk -F: '{print $2}' | sed 's/ //g'`
name_split="$run_name.genfastq.$rnds"
name_cjobs="$run_name.create_jobs.$rnds" 
name_cmd="$run_name.submit_jobs.$rnds"

to_lsf $name_split "$prod_pipe/bin/bfast.split.reads.rb $config_file"
to_lsf $name_cjobs "$prod_pipe/bin/bfast.lsf.submit.rb $config_file" $name_split
to_lsf $name_cmd "./cluster_JOBS.sh $config_file" $name_cjobs
