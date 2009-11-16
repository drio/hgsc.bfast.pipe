#!/usr/bin/env ruby19

%w(yaml fileutils).each { |dep| require dep }

# Deals with the creation and deps of LSF jobs
class LSFDealer
  def initialize(seed, queue, output_file="cluster_JOBS.sh")
    @queue = queue
    @output_file = output_file
    @contents = [] # contents of the job_list script
    @contents << "#!/bin/bash\n\n"
		@n_jobs = 1
    @log_dir = "./lsf_logs"
    # We'll use this for the job name
		@seed = seed
    FileUtils.mkdir_p @log_dir
  end

  # Add a regular job
  def add_job(job_root, cmd, deps=nil)
    job_name = @seed + "." + job_root + "." +  @n_jobs.to_s + rand(10000).to_s
    wait_for = deps.nil? ? "" : (find_deps deps)
    @contents << "bsub -o #{@log_dir}/#{job_name}.out \\"
    @contents << "-e #{@log_dir}/#{job_name}.err -J #{job_name} \\"
    @contents << "-q #{@queue} #{wait_for} \"#{cmd}\""
		@n_jobs += 1
		job_name
  end

	def blank
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

# Reads the bfast experiment config
class Config
	def initialize(config)
	 %w(input global match local post).each {|r| set config, r }
	end

	def set(config, r)
		config["#{r}_options"].each do |key, value| 
			instance_variable_set("@#{r}_#{key}", value)
		end
	end
end

#
# Main
# Load config
config = Config.new( YAML::load(DATA) )

# Prepare LSF script
lsf = LSFDealer.new("run-name", "normal")

# split job
dep1 = []
dep1 << lsf.add_job("job1", "xxx1")
lsf.blank

# Per each fastq read file:
# match/localalign/postprocess/bam/sort
dep2 = []
dep2 << lsf.add_job("job2", "xxx2", dep1)
dep2 << lsf.add_job("job3", "xxx3", dep1)
lsf.blank

# Merge all the previous BAMs
lsf.add_job("job4", "xxx3", dep2)
lsf.create_file


__END__
input_options:
 run_name: run_name
 f3: reads_f3
 r3: reads_r3
 qf3: quals_f3
 qr3: quals_r3
global_options:
 bfast_bin: /usr/local/bfast/bin
 samtools_in: /usr/local/samtools
 space: CS
 fasta_file_name: /genomes/hg18/hg18.fasta
 timing: ON
 run_directory: /rundir
 reads_directory: /reads_dir
 output_directory: /output_dir
 tmp_directory: /tmp_dir
 output_id: output_id_like_test_drio
 splits: 10
match_options:
 threads: 8
local_options:
 threads: 8
post_options:
 algorithm: 4
