# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

%w(yaml fileutils singleton digest/md5).each { |dep| require dep }

# Deals with the creation and deps of LSF jobs
class LSFDealer
  def initialize(seed, queue, wrap_cmd_scpt, t_dir, n_splits=1, output_file=nil)
    @queue = queue
    @output_file ||= "cluster_JOBS.sh"
    @contents = [] # contents of the job_list script
    @contents << "#!/bin/bash\n\n"
    @n_jobs = 1
    @log_dir = "./lsf_logs"
    @track_dir = t_dir
    # We'll use this for the job name
    @seed = seed
    # Check if the wrap_cmd_scpt is there
    File.exist?(wrap_cmd_scpt) or raise "wrap_cmd_scpt not found."
    @wrap_cmd_scpt = wrap_cmd_scpt
    FileUtils.mkdir_p @log_dir
    (1..n_splits).each {|s| FileUtils.mkdir_p "./#{@log_dir}/#{s}" }
  end

  # Add a regular job
  def add_job(job_root, cmd, split_n=1, resources=nil, deps=nil)
    r_split  = split_n == "" ? "" : ".split"
    job_name = @seed + "." + job_root + r_split + split_n.to_s + "." + 
               random_string
               # + "." + @n_jobs.to_s + rand(100).to_s
    wait_for = deps.nil? ? "" : (find_deps deps)

    # If jobs already run successfully, we don't want to run it
    if File.exists?(@track_dir + "/" + job_name + ".ok")
      puts "skipping #{job_name}"
      return ""
    else
      puts "Adding #{job_name}"
      @contents << "echo 'Submitting: #{job_name}'"
      @contents << "bsub -o #{@log_dir}/#{split_n}/#{job_name}.out \\"
      @contents << "-e #{@log_dir}/#{split_n}/#{job_name}.err \\"
      @contents << "-J #{job_name} \\"
      @contents << "-q #{@queue} #{wait_for} \\"
      @contents << "-R '#{resources}' \\" if resources
      @contents << wrap_cmd(job_name, cmd)
      @contents << ""
      @n_jobs += 1
      job_name
    end
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

  def random_string
    Digest::MD5.hexdigest(rand(10000000000).to_s)[1..5]
  end

  # Wrap the cmd around script.run.process
  def wrap_cmd(job_name, cmd)
    " \"#{@wrap_cmd_scpt} #{job_name} '#{cmd}'\" "
  end

  # Gets the array of deps you want to associate to a job
  def find_deps(s_deps)
    list_of_deps = "-w '"
    s_deps.each do |d|
      list_of_deps << "done(#{d})" + " && " if d != ''
    end
    list_of_deps.gsub!(/ && $/,"'")
  end
end
