%w(yaml fileutils singleton).each { |dep| require dep }

# Miscelanea
module Misc
	VERSION = "0.1a"
end

# Deals with the User Interface
class UInterface
	include Singleton

	# Loads config as argument or uses DATA
	def load_config(arguments, tool)
		@tool = tool
		(arguments.size == 1) ? File.new(arguments[0]) : error("Config not found")
	end

	def error(msg)
		puts "Error: " + msg; puts "" 
    puts usage
		exit 1
	end

	private

	def usage
		template = %Q{
		bfast hgsc pipeline
    VERSION: xversionx

    Usage: xcmdx <config_file>
		}

    template.gsub!(/xversionx/, Misc::VERSION)
    template.gsub!(/xcmdx/    , File.basename(@tool))
    template.gsub!(/^\s+/     , '')
	end
end

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
		"#{samtools} sort #{bam_file} #{bam_file_sorted} "
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

# Reads the bfast experiment config
class Config
	def initialize(config)
	 %w(input global match local post).each {|r| set config, r }
	end

	# Set all the config entries as methods for this class
	def set(config, r)
		config["#{r}_options"].each do |key, value|
			instance_variable_set("@#{r}_#{key}", value)
		end
	end

	# If method missing, use the method name and return the value of the instance
  # variable with that name
	def method_missing(m, *args, &block)
		eval "@#{m}"
	end
end

# Deals with the creation and deps of LSF jobs
class LSFDealer
  def initialize(seed, queue, n_splits, output_file="cluster_JOBS.sh")
    @queue = queue
    @output_file = output_file
    @contents = [] # contents of the job_list script
    @contents << "#!/bin/bash\n\n"
		@n_jobs = 1
    @log_dir = "./lsf_logs"
    # We'll use this for the job name
		@seed = seed
    FileUtils.mkdir_p @log_dir
		(1..n_splits).each {|s| FileUtils.mkdir_p "./#{@log_dir}/#{s}" }
  end

  # Add a regular job
  def add_job(job_root, cmd, split_n, resources=nil, deps=nil)
    job_name = @seed + "." + job_root + ".split" + split_n.to_s + "." +
               @n_jobs.to_s + rand(100).to_s
    wait_for = deps.nil? ? "" : (find_deps deps)
    @contents << "bsub -o #{@log_dir}/#{split_n}/#{job_name}.out \\"
    @contents << "-e #{@log_dir}/#{split_n}/#{job_name}.err \\"
    @contents << "-J #{job_name} \\"
    @contents << "-q #{@queue} #{wait_for} \\"
    @contents << "-R '#{resources}' \\" if resources
    @contents << "\"#{cmd}\""
    @contents << ""
		@n_jobs += 1
		job_name
  end

	def blank(msg=nil)
    @contents << "# " + msg if msg
    @contents << ""
	end

  # Dump the contents to a file
  def create_file
    File.open("./#{@output_file}", "w") { |f| f.puts(@contents) }
    FileUtils.chmod 0755, "./#{@output_file}"
  end

	private

  # Gets the array of deps you want to associate to a job
  def find_deps(s_deps)
    list_of_deps = "-w '"
    s_deps.each { |d| list_of_deps << "done(#{d})" + " && " }
    list_of_deps.gsub!(/ && $/,"'")
  end
end

