#!/bin/bash
#
# Runs the ruby tool to generate csv output of our 
# SOLiD data. It checks the exit code and sends and
# email.
#
help()
{
	echo "$1"
  cat <<EOF
Usage:
  $0 cmd file_to_dump_csv log_file mail_to
Example:
  $0 "./bin/list_system_data.rb /stornext/snfs4/next-gen/solid/analysis" /tmp/foo.csv /tmp/log.txt "deiros@bcm.edu dc12@bcm.edu"
EOF
	exit 1
}

send_email()
{
  error_ok=$1
  subject="csv dump tool heartbeat: {$error_ok}"

  if [ $error_ok == "OK" ] 
  then
    n_seas=$[$(cat $dump_to | wc -l) - 1]
  else
    n_seas="-"
  fi

  (
  cat <<EOF
Heart Beat for the CSV dumping of our SOLiD data
------------------------------------------------

Started : $started
Finished: $finished

# of seas processed: $n_seas

csv file: $dump_to
log file: $log_file
EOF
  ) | mail -s "$subject" $mail_to
}

######################################################
# Main
######################################################
[ ".$1" == "." ] && help "Error param: cmd"
[ ".$2" == "." ] && help "Error param: dump_to file"
[ ".$3" == "." ] && help "Error param: log file"
[ ".$4" == "." ] && help "Error param: mail_to"

cmd="$1"
dump_to="$2"
log_file="$3"
mail_to="$4"

started=`date`
$cmd > $dump_to 2> $log_file &
wait
exit_code=$?
finished=`date`

if [ $exit_code -eq 0 ]
then
  send_email "OK"
	exit 0
else
  send_email "ERROR"
	exit 1
fi
