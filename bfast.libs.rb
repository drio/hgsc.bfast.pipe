%w(yaml fileutils).each { |dep| require dep }

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

