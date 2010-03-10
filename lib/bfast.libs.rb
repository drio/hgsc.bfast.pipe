# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

%w(yaml fileutils singleton digest/md5).each { |dep| require dep }

# Miscelanea
module Misc
  VERSION  = "0.1"
  RG_FNAME = "./rg.txt"

  # Create a file with the @RG line that bfast will use to build the @RG type
  # in the header of the bam
  #rg_id: 0                                 ; Read group (BAM spec) 
  #rg_pl: SOLiD                             ; platform 
  #rg_pu: SOLiD0097_20081114_1_Hs_1011_MP_F3; Our uniq identifier RUN + LIMS ID
  #rg_lb: HS_1011_MP                        ; Library 
  #rg_ds: rl=25                             ; Description rl=25 OR rl=50_25
  #rg_dt: 2010-01-22T18:20:29-0600          ; date, ISO 8601
  #rg_sm: CMT-001                           ; sample / pool name
  #rg_cn: Baylor                            ; center
  def Misc.create_rg_file(config)
    rg_data = "@RG\t"
    %w{ID PL PU LB DS DT SM CN}.each do |t|
      t_value = eval("config.post_rg_" + t.downcase).to_s
      puts "#{t}:-#{t_value}--" 
      raise "Couldn't find #{t} tag in config file." unless t_value
      rg_data << "#{t}:" + t_value + "\t"
    end

    File.open(RG_FNAME, "w") {|f| f << rg_data}
  end

  # Figure out how to flag bfast to process input data
  def Misc.input_compress(config)
    case config.global_compress_input
      when "bzip2"
        "-j "
      when "gzip"
        "-z "
      else  
        " "
    end
  end

  # Figure out how to flag bfast to process output data
  def Misc.output_compress(config)
    case config.global_compress_splits
      when "bzip2"
        "-J "
      when "gzip"
        "-Z "
      else  
        " "
    end
  end

  # Generate the wildcard to find files based on the user preferences
  def Misc.wild(root, config)
    case config.global_compress_input
      when "bzip2"
        "*.#{root}.bz2 "
      when "gzip"
        "*.#{root}.gz "
      else  
        "*.#{root}"
    end
  end
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
    main_bin('match') + core_cmd + Misc::input_compress(@config) +
    " #{tmp_arg} -r #{fastq} > " + match_file
  end
  
  # bfast localalign -A 1 -t -n 8 -f $ref
  # -m bfast.matches.file.$root.bmf # bfast.aligned.file.$root.baf
  def local
    main_bin('localalign') + core_cmd + " -m #{match_file} > #{local_file}"
  end
  # -a 3 -O 3 > bfast.reported.file.$root.sam
  def post
    main_bin('postprocess') + " -f #{ref} " + " -i #{local_file} " +
    "-a 3 -O 3 -r #{Misc::RG_FNAME} > #{sam_file}"
    #"-n 8 -a 3 -O 3 > #{sam_file}"
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
    #"#{samtools} sort #{bam_file} #{bam_file_sorted} "
    "#{@config.global_java_vm} -jar -Xmx#{java_vm_mem_sort} " +
    "#{picardjars}/SortSam.jar " +
    "TMP_DIR=#{@config.global_tmp_dir} " +
    "INPUT=#{merged_bam} " +
    "OUTPUT=#{bam_file_sorted} " +
    "SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT"
  end

  def dups
    "#{@config.global_java_vm} -jar -Xmx#{java_vm_mem_dups} " +
    "#{picardjars}/MarkDuplicates.jar " +
    "TMP_DIR=#{@config.global_tmp_dir} " +
    "INPUT=#{bam_file_sorted} " +
    "OUTPUT=#{bam_file_sorted_dups} " +
    "METRICS_FILE='./metric_file.picard' " +
    "VERBOSITY=ERROR " +
    "VALIDATION_STRINGENCY=SILENT "
  end

  # samtools index bfast.reported.file.$root.sorted.bam
  def index2
    "#{samtools} index #{bam_file_sorted} "
  end

  # samtools merge out.bam in1.bam in2.bam in3.bam
  def final_merge
    #"#{samtools} merge #{merged_bam} #{list_bams_to_merge}"o
    "#{@config.global_java_vm} -jar -Xmx#{java_vm_mem_dups} " +
    "#{picardjars}/MergeSamFiles.jar " +
    "#{list_bams_to_merge} " +
    "TMP_DIR=#{@config.global_tmp_dir} " +
    "OUTPUT=#{merged_bam} " +
    "VALIDATION_STRINGENCY=SILENT "
  end

  # Regenerate the header so we have more useful information on it
  def gen_header
    cmd = "#{@config.global_java_vm} -jar #{@config.header_regen_jar} "
    cmd << "type=" + @config.header_sq_type + " "
    #%w{ID PL PU LB DS DT SM CN}.each do |t|
    # t_value = eval("@config.post_rg_" + t.downcase).to_s
    # cmd << "#{t}=" + t_value + " "
    #end
    cmd << "I=#{bam_file_sorted_dups} "
    cmd << "O=#{bam_file_sorted_dups_fix_header} "
  end

  # computes the stats
  def stats_frag
    stats_core + " 3 solid > marked.stats.txt"
  end

  def stats_f3
    stats_core + " 1 solid > marked.stats.F3.txt"
  end

  def stats_r3
    stats_core + " 2 solid > marked.stats.R3.txt"
  end

  def capture_stats
    "#{@config.global_java_vm} -cp #{@config.capture_j_classpath}" +
    " -Xmx6000M CaptureStatsBAM4 -o #{root_name} -t " +
    "#{@config.capture_chip_design} " +
    "-i #{bam_file_sorted_dups} -w -d"
  end

  def bam_reads_validator
    read_val_core + " #{bam_file_sorted_dups} #{@config.global_run_dir}" +
    " ./output/#{root_name}.read_val_log.txt"
  end

  def email_success
    email_to  = @config.success_email_to
    cat_files = if @config.global_input_MP == 1
      "marked.stats.F3.txt marked.stats.R3.txt"
    else 
      "marked.stats.txt"
    end
    cmd = "cat #{cat_files} | "
    cmd << 'mail -s \"[OK] BFast analysis completed: ' + root_name + '\"'
    cmd << " #{email_to}"
  end

  private
  
  def stats_core
    "#{@config.global_java_vm} -jar #{@config.stats_s_jar} #{bam_file_sorted_dups} "
  end

  def read_val_core
    "#{@config.global_java_vm} -jar #{@config.stats_bam_reads_val_jar} "
  end

  def java_vm_mem_dups; "#{@config.dups_java_vm_mem}"; end

  def java_vm_mem_sort; "#{@config.sort_java_vm_mem}"; end

  def merged_bam; "./output/#{root_name}.merged.bam"; end

  def list_bams_to_merge
    t = "output/SS/bfast.reported.file.#{root_name}.SS.bam"
    (1..@n_splits).inject("") {|list, s| list << "INPUT=" + t.gsub(/SS/, "split#{s}") + " " }
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
     "-A 1 -t -n #{@config.global_threads} -f #{ref} "
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
    "#{@config.global_samtools_bin}/samtools"
  end

  def picardjars
    "#{@config.global_picardjars}/"
  end

  def bam_file
    split_dir + "/bfast.reported.file.#{root_name}.#{@curr_split}.bam"
  end

  def bam_file_sorted
    "./output/#{root_name}.sorted.bam"
  end

  def bam_file_sorted_dups
    "./output/#{root_name}.sorted.dups.bam"
  end

  def bam_file_sorted_dups_fix_header
    "./output/#{root_name}.sorted.dups.with.header.bam"
  end

  def set_current_split(fastq_file)
    @curr_split = "split" + fastq_file.split(".")[-2]
  end
end

# Reads the bfast experiment config yaml and prepares config
class Config
  def initialize(config)
    %w(input global match local post tobam
       sort dups final header stats capture success).each do |r|
      set config, r
    end
  end

  # Set all the config entries as methods for this class
  def set(config, r)
    puts "Loading config option: #{r}_options"
    config["#{r}_options"].each do |key, value|
      puts "loading suboption: #{r}_#{key}"
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

