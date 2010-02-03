#!/usr/bin/env ruby19

$: << File.join(File.dirname(File.dirname($0)), "lib")
%w(bfast.libs).each { |dep| require dep }

ui = UInterface.instance
config = Config.new( YAML::load(ui.load_config(ARGV, $0)) )
cmd = "#{config.global_bfast_bin}/../scripts/solid2fastq " +
      "-n #{config.global_reads_per_file} " +
      "-o #{config.input_run_name} " +
      Misc::input_compress(config) +
      Misc::output_compress(config) +
      "#{config.global_run_dir}/" + Misc::wild("csfasta", config) + " " +
      "#{config.global_run_dir}/" + Misc::wild("qual", config)
puts cmd
Dir.chdir(config.global_reads_dir) { `#{cmd}` }

__END__
input_options:
 run_name: run_small_test
 f3: reads_f3
 r3: reads_r3
 qf3: quals_f3
 qr3: quals_r3
global_options:
 bfast_bin: /stornext/snfs1/next-gen/drio-scratch/bfast_related/bfast/bfast
 samtools_in: /stornext/snfs1/next-gen/software/samtools-0.1.6
 space: CS
 fasta_file_name: /genomes/hg18/hg18.fasta
 timing: ON
 logs_dir: /stornext/snfs3/drio_scratch/small_test/logs
 run_dir: /stornext/snfs3/drio_scratch/small_test/input
 reads_dir: /stornext/snfs3/drio_scratch/small_test/reads
 output_dir: /stornext/snfs3/drio_scratch/small_test/output
 tmp_dir: /tmp_dir
 output_id: output_id_like_test_drio
 reads_per_file: 25000
match_options:
 threads: 8
local_options:
 threads: 8
post_options:
 algorithm: 4
