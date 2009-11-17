#!/usr/bin/env ruby19

$:.unshift File.join(File.dirname(__FILE__))
%w(bfast.libs).each { |dep| require dep }

# Encapsulates the creation of bfast commands
class BfastCmd
	def initialize(config)
		@config = config
	end

	# bfast match -A 1 -t -n 8 -f $ref -r $fastq > bfast.matches.file.$root.bmf
	def match(fastq)
		main_bin('match') + core_cmd + " -r #{fastq} > " + match_file
	end
	
	# bfast localalign -A 1 -t -n 8 -f $ref
  # -m bfast.matches.file.$root.bmf # bfast.aligned.file.$root.baf
	def local
		main_bin('local') + core_cmd + " -m #{match_file} > #{local_file}"
	end

	# bfast postprocess -f $ref -i bfast.aligned.file.$root.baf 
  # -a 3 -O 3 > bfast.reported.file.$root.sam
	def post
		main_bin('postprocess') + " -f #{ref} " + " -i #{local_file} " +
		"-a 3 -O 3 > #{sam_file}"
	end

	# samtools view -bt $ref 
  # -o bfast.reported.file.$root.bam bfast.reported.file.$root.sam
	def tosam
		"#{samtools} view -bt #{ref} -o #{bam_file} #{sam_file}"
	end

	# samtools index bfast.reported.file.$root.bam
	def index1
		"#{samtools} index #{bam_file} "
	end

  # samtools sort ./bfast.reported.file.$root.bam 
  # ./bfast.reported.file.$root.sorted
	def sort
		"#{samtools} sort #{bam_file} -o #{bam_file_sorted} "
	end

  # samtools index bfast.reported.file.$root.sorted.bam
	def index2
		"#{samtools} index #{bam_file_sorted} "
	end

	# samtools merge out.bam in1.bam in2.bam in3.bam
	def final_merge
		"#{samtools} merge #{root_name}.merged.bam *.sorted.bam"
	end	

	private
	
	def main_bin(sub_cmd)
		"#{@config.global_bfast_bin}/bfast #{sub_cmd} "
	end

	def core_cmd
		"-A 1 -t -n 8 -f #{ref}"
	end

	def match_file
		"bfast.matches.file.#{root_name}.bmf"
	end

	def local_file
		"bfast.matches.file.#{root_name}.baf"
	end

	def ref
		@config.global_fasta_file_name
	end

	def sam_file
		"bfast.reported.file.#{root_name}.sam"
	end

	def root_name
		@config.input_run_name
	end

	def samtools
		"#{@config.global_samtools_bin}/samtools "
	end

	def bam_file
		"bfast.reported.file.#{root_name}.bam"
	end

	def sam_file
		"bfast.reported.file.#{root_name}.sam"
	end

	def bam_file_sorted
		"bfast.reported.file.#{root_name}.sorted.bam"
	end
end

# Load config
config = Config.new( YAML::load(DATA) )
# Prepare LSF 
lsf    = LSFDealer.new(config.input_run_name, config.global_lsf_queue)
# Prepare bfast cmd generation
cmds   = BfastCmd.new(config)

# Get list of read splits/files
splits = Dir[config.global_reads_dir + "/*.fastq"]

# Per each split, create the basic bfast workflow with deps
final_deps = []
splits.each do |s|
	dep = lsf.add_job("match" , cmds.match(s))
	dep = lsf.add_job("local" , cmds.local , [dep])
	dep = lsf.add_job("postp" , cmds.post  , [dep])
	dep = lsf.add_job("tosam" , cmds.tosam , [dep])
	dep = lsf.add_job("index1", cmds.index1, [dep])
	dep = lsf.add_job("sort"  , cmds.sort  , [dep])
	dep = lsf.add_job("index2", cmds.index2, [dep])
	lsf.blank

	final_deps << dep
end

# when all the previous jobs are completed, we can merge all the bams
lsf.add_job("final_merge", cmds.final_merge, final_deps)

lsf.create_file

__END__
input_options:
 run_name: run_small_test
 f3: reads_f3
 r3: reads_r3
 qf3: quals_f3
 qr3: quals_r3
global_options:
 lsf_queue: normal
 bfast_bin: /stornext/snfs1/next-gen/drio-scratch/bfast_related/bfast
 samtools_bin: /stornext/snfs1/next-gen/software/abi/bioscope/corona/bin
 space: CS
 fasta_file_name: /stornext/snfs1/next-gen/solid/cmaps/h/human_g1k_v37/human_g1k_v37.fa
 timing: ON
 logs_dir: /stornext/snfs3/drio_scratch/small_test/logs
 run_dir: /stornext/snfs3/drio_scratch/small_test/input
 reads_dir: /stornext/snfs3/drio_scratch/small_test/reads
 output_dir: /stornext/snfs3/drio_scratch/small_test/output
 tmp_dir: /space1/tmp/tmp.bfast.solid.pipe
 reads_per_file: 100000
match_options:
 threads: 8
local_options:
 threads: 8
post_options:
 algorithm: 4
