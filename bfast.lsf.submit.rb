#!/usr/bin/env ruby19

$:.unshift File.join(File.dirname(__FILE__))
%w(bfast.libs).each { |dep| require dep }

# Encapsulates the creation of bfast commands
class BfastCmd
	def initialize(config, read_files)
		@config   = config
		@n_splits = read_files.size
		create_output_dirs
	end

	def create_output_dirs
		(1..@n_splits).each {|s| FileUtils.mkdir_p "./output/split#{s}" }
	end

	# bfast match -A 1 -t -n 8 -f $ref -r $fastq > bfast.matches.file.$root.bmf
	def match(fastq)
		set_current_split(fastq)
		main_bin('match') + core_cmd + " #{tmp_arg} -r #{fastq} > " + match_file
	end
	
	# bfast localalign -A 1 -t -n 8 -f $ref
  # -m bfast.matches.file.$root.bmf # bfast.aligned.file.$root.baf
	def local
		main_bin('localalign') + core_cmd + " -m #{match_file} > #{local_file}"
	end

	# bfast postprocess -f $ref -i bfast.aligned.file.$root.baf 
  # -a 3 -O 3 > bfast.reported.file.$root.sam
	def post
		main_bin('postprocess') + " -f #{ref} " + " -i #{local_file} " +
		"-a 3 -O 3 > #{sam_file}"
	end

	# samtools view -bt $ref 
  # -o bfast.reported.file.$root.bam bfast.reported.file.$root.sam
	def tobam
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
		"#{samtools} merge ./output/#{root_name}.merged.bam #{list_bams_to_merge}"
	end

	private

	def list_bams_to_merge
		t = "output/SS/bfast.reported.file.#{root_name}.SS.sorted.bam"
		(1..@n_splits).inject("") {|list, s| list << t.gsub(/SS/, "split#{s}") + " " }
	end

	def split_dir
		"./output/#{@curr_split}"
	end
	
	def main_bin(sub_cmd)
		"#{@config.global_bfast_bin}/bfast #{sub_cmd} "
	end

	def tmp_arg
		"-T #{@config.global_tmp_dir}"
	end

	def core_cmd
		 "-A 1 -t -n 8 -f #{ref}"
	end

	def match_file
		split_dir + "/bfast.matches.file.#{root_name}.#{@curr_split}.bmf"
	end

	def local_file
		split_dir + "/bfast.matches.file.#{root_name}.#{@curr_split}.baf"
	end

	def ref
		@config.global_fasta_file_name
	end

	def sam_file
		split_dir + "/bfast.reported.file.#{root_name}.#{@curr_split}.sam"
	end

	def root_name
		@config.input_run_name
	end

	def samtools
		"#{@config.global_samtools_bin}/samtools "
	end

	def bam_file
		split_dir + "/bfast.reported.file.#{root_name}.#{@curr_split}.bam"
	end

	def bam_file_sorted
		split_dir + "/bfast.reported.file.#{root_name}.#{@curr_split}.sorted"
	end

	def set_current_split(fastq_file)
		@curr_split = "split" + fastq_file.split(".")[-2]
	end
end

# Load config
config = Config.new( YAML::load(DATA) )
# Prepare LSF 
lsf    = LSFDealer.new(config.input_run_name, config.global_lsf_queue)
# Get list of read splits/files
splits = Dir[config.global_reads_dir + "/*.fastq"]
# Prepare bfast cmd generation
cmds   = BfastCmd.new(config, splits)

# Per each split, create the basic bfast workflow with deps
one_machine = "rusage[mem=29000]span[hosts=1]"
reg_job     = "rusage[mem=4000]"
final_deps = []
splits.each do |s|
	dep = lsf.add_job("match" , cmds.match(s), one_machine)
	dep = lsf.add_job("local" , cmds.local , one_machine, [dep])
	dep = lsf.add_job("postp" , cmds.post  , reg_job    , [dep])
	dep = lsf.add_job("tobam" , cmds.tobam , reg_job    , [dep])
	dep = lsf.add_job("index1", cmds.index1, reg_job    , [dep])
	dep = lsf.add_job("sort"  , cmds.sort  , reg_job    , [dep])
	dep = lsf.add_job("index2", cmds.index2, reg_job    , [dep])
	lsf.blank "----------------"

	final_deps << dep
end

# when all the previous jobs are completed, we can merge all the bams
lsf.add_job("final_merge", cmds.final_merge, reg_job, final_deps)

lsf.create_file

__END__
input_options:
 run_name: run_small_test
 f3: reads_f3
 r3: reads_r3
 qf3: quals_f3
 qr3: quals_r3
global_options:
 lsf_queue: test
 bfast_bin: /stornext/snfs1/next-gen/drio-scratch/bfast_related/bfast/bfast
 samtools_bin: /stornext/snfs1/next-gen/software/samtools-0.1.6
 space: CS
 fasta_file_name: /stornext/snfs3/drio_scratch/bf.indexes/small/test.fasta
 timing: ON
 logs_dir: /stornext/snfs3/drio_scratch/small_test/logs
 run_dir: /stornext/snfs3/drio_scratch/small_test/input
 reads_dir: /stornext/snfs3/drio_scratch/small_test/reads
 output_dir: /stornext/snfs3/drio_scratch/small_test/output
 tmp_dir: /space1/tmp/
 reads_per_file: 100000
match_options:
 threads: 8
local_options:
 threads: 8
post_options:
 algorithm: 4
