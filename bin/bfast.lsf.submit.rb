#!/usr/bin/env ruby19

$: << File.join(File.dirname(File.dirname($0)), "lib")
require 'load_libs'

# Load config
ui = UInterface.instance
config = Config.new( YAML::load(ui.load_config(ARGV, $0)) )

# Find path to the wrapper script to run bfast cmds and check 
# for successful execution
cmd_wrapper_scpt = File.dirname(__FILE__) + "/script.run.process.sh"

# Get list of read splits/files
#splits = Dir[config.global_reads_dir + "/*.fastq"]
path_to_fastqs = config.global_reads_dir + "/" + 
                 Misc::wild("fastq", config).gsub(/ /,'')
splits = Dir[path_to_fastqs.chomp]
puts "Splits wild = -#{path_to_fastqs.chomp}-"
puts "Splits found = #{splits.size}"

# Prepare LSF 
lsf = LSFDealer.new(config.input_run_name,
                    config.global_lsf_queue,
	  								cmd_wrapper_scpt,
	  								config.global_trackdir,
                    splits.size)
# Prepare bfast cmd generation
cmds = BfastCmd.new(config, splits)

# Create the rg.txt file (@RG tag)
Misc::create_rg_file(config)

# Per each split, create the basic bfast workflow with deps
#reg_job     = "rusage[mem=4000]"
#one_machine = "rusage[mem=28000]span[hosts=1]"
re_match  = config.match_lsf_resources
re_local  = config.local_lsf_resources
re_post   = config.post_lsf_resources
re_tobam  = config.tobam_lsf_resources
re_sort   = config.sort_lsf_resources
re_dups   = config.dups_lsf_resources
re_final  = config.final_lsf_resources
re_header = config.header_lsf_resources
re_stats  = config.stats_lsf_resources
re_cap    = config.capture_lsf_resources

final_deps = []
splits.each do |s|
	sn  = s.match(/\.(\d+)\./)[1]
	puts "Jobs for split: #{s.split} - #{sn}"
	dep = lsf.add_job("match" , cmds.match(s), sn, re_match)
	dep = lsf.add_job("local" , cmds.local   , sn, re_local , [dep])
	dep = lsf.add_job("postp" , cmds.post    , sn, re_post  , [dep])
	dep = lsf.add_job("tobam" , cmds.tobam   , sn, re_tobam , [dep])
	lsf.blank "----------------"

	final_deps << dep
end

# when all the previous jobs are completed, we can merge all the bams
dep = lsf.add_job("merge" , cmds.final_merge, "", re_final, final_deps)

# Sort and mark dups in the final BAM
dep = lsf.add_job("sort", cmds.sort, "", re_sort, [dep])
dep = lsf.add_job("dups", cmds.dups, "", re_dups, [dep])
# Regenerate the header so we have better and more clear @SQ entries
dep = lsf.add_job("regen_bam_header", cmds.gen_header, "", re_header, [dep])

# Run stats
s_deps = []
if config.global_input_MP == 0
  s_deps << lsf.add_job("stats", cmds.stats_frag, "", re_stats, [dep])
else
  s_deps << lsf.add_job("stats_F3", cmds.stats_f3, "", re_stats, [dep])
  s_deps << lsf.add_job("stats_R3", cmds.stats_r3, "", re_stats, [dep])
end

# Run BAM Reads Validation
dep = lsf.add_job("bam_reads_val", cmds.bam_reads_validator, "", re_stats, [dep])

# Run Capture Stats
caps_dir = config.capture_stats_dir
if config.global_input_CAP == 1
  Dir.mkdir(caps_dir)
  dep = lsf.add_job("capture_stats", cmds.capture_stats, "", re_cap, s_deps)
end

# Email if the analysis went well
email_deps = config.global_input_CAP == 0 ? s_deps : [dep]
lsf.add_job("email_success", cmds.email_success, "", nil, email_deps)

lsf.create_file

__END__
input_options:
 run_name: run_small_test
 f3: reads_f3
 r3: reads_r3
 qf3: quals_f3
 qr3: quals_r3
global_options:
 lsf_queue: hptest
 threads: 7
 bfast_bin: /stornext/snfs1/next-gen/drio-scratch/bfast_related/bfast/bfast
 samtools_bin: /stornext/snfs1/next-gen/software/samtools-0.1.6
 space: CS
 fasta_file_name: /stornext/snfs3/drio_scratch/bf.indexes/small/test.fasta
 timing: ON
 logs_dir: ./lsf_logs
 run_dir: /BA_EVAL_01/drio_scratch/bfast/small_test/input
 reads_dir: /BA_EVAL_01/drio_scratch/bfast/small_test/reads
 output_dir: /BA_EVAL_01/drio_scratch/bfast/small_test/output
 tmp_dir: /space1/tmp/
 output_id: output_id_like_test_drio
 reads_per_file: 50000
match_options:
 threads: 8
 lsf_resources: "rusage[mem=2800]span[hosts=1]"
local_options:
 threads: 8
 lsf_resources: "rusage[mem=2800]span[hosts=1]"
post_options:
 algorithm: 4
 lsf_resources: "rusage[mem=400]"
 rg_id: 0
 rg_pl: solid
 rg_pu: SOLiD0097_20081114_1_Hs_1011_MP_F3
 rg_lb: HS_1011_MP
 rg_ds: rl=25
 rg_dt: 2010-01-22T18:20:29-0600
 rg_sm: CMT-001
 rg_cn: Baylor
tobam_options:
 lsf_resources: "rusage[mem=400]"
final_options:
 lsf_resources: "rusage[mem=400]"
