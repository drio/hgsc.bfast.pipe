#!/usr/bin/perl -w

##!/data/pipeline/code/production/bin/x86_64-linux/perl -w

# EXAMPLE to Run: "perl ./getLibraryInfo.pl 0000001003220613-1

#INPUT:
# barcode
# - barcode of the Spot, which is the slide barcode plus the addition of the spot number (ie, then the spot barcode is 0000001003220613-2).
# 
# Here is the format for the output:
# Library_id;patient_id;sample_id;project_name
# 
# Here is what I have filled out for the rg fields:
#  rg_id: 0 (always 0)
#  rg_pl: SOLiD (always SOLiD)
#  rg_pu: 0301_20091229_2_SL_AWG_LSBC_3204_T_000sA_01003244678_1 (full run
# name)
#  rg_lb: AWG_LSBC.3204.T_000sA (lib_id)
#  rg_ds: rl=50 (length)
#  rg_dt: 2010-02-07T00:43:20 (currently dt is when I have started the
# analysis.)
#  rg_sm: 3204.T (sample_id)
#  rg_cn: Baylor (always Baylor)


use strict;
use LWP;

my $ncbiURL ="http://10.10.53.190:8080/solidlims/getLibraryInfo.jsp?";
my $paraStr = "barcode=" . $ARGV[0];
my $aTag = "";
if ( defined($ARGV[1]) ) {
     $aTag ="tag=" . $ARGV[1];
}

$ncbiURL="$ncbiURL$paraStr$aTag";
#print "$ncbiURL\n";

my $ua = LWP::UserAgent->new;
my $response=$ua->get($ncbiURL);


if(not $response->is_success ) {print "Error: Cannot connect\n"; exit(0);}

my $textStr= $response->content;
$textStr=~/^\s*(.+)\s*$/;
$textStr=$1;
print "$textStr\n";
