require 'tempfile'
# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

%w(yaml fileutils singleton digest/md5).each { |dep| require dep }

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

  # Flag the DB so we know the analysis has started
  def fdb_completed
    update_db_tool          +
    '-a update '            +
    '-n "#{root_name}" '    +
    "--key=completed "      +
    '--value="#{Time.now}"'
  end

  # bfast match -A 1 -t -n 8 -f $ref -r $fastq > bfast.matches.file.$root.bmf
  def match(fastq)
    set_current_split(fastq)
    main_bin('match') + core_cmd + Misc::input_compress(@config) +
    " #{tmp_arg} -r #{fastq} > " + match_file
  end

  # Use bwaaln instead of the match step (experimental)
  # TODO: Eventually, bfast2_bin will change to main_bin 
  def bwaaln(fastq)
    set_current_split(fastq)
    bfast2_bin('bwaaln') + "-c -t#{@config.global_threads} #{@config.bwaaln_prefix} " + 
    "#{fastq} > " + match_file
  end
  
  # bfast localalign -A 1 -t -n 8 -f $ref
  # -m bfast.matches.file.$root.bmf # bfast.aligned.file.$root.baf
  def local
    main_bin('localalign') + core_cmd + " -m #{match_file} > #{local_file}"
  end

  # Since we match with bwaaln, we have to enable -U:
  # -U    Do not use mask constraints from the match step
  # TODO Eventually, bfast2_bin will change to main_bin 
  def local_u
    bfast2_bin('localalign') + core_cmd + " -U -m #{match_file} > #{local_file}"
  end

  # -a 3 -O 3 > bfast.reported.file.$root.sam
  def post
    main_bin('postprocess') + " -f #{ref} " + " -i #{local_file} " +
    "-a 3 -O 1 > #{sam_file}"

    ## If you want RG tag
    ##"-a 3 -O 3 -r #{Misc::RG_FNAME} > #{sam_file}"

    #"-n 8 -a 3 -O 3 > #{sam_file}"
  end

  # samtools view -bt $ref 
  # -o bfast.reported.file.$root.bam bfast.reported.file.$root.sam
  def tobam
    #"#{samtools} view -bt #{ref} -o #{bam_file} #{sam_file}"
    "#{@config.global_java_vm} -jar -Xmx#{java_vm_mem_to_bam} " +
    "#{picardjars}/SamFormatConverter.jar " +
    "TMP_DIR=#{@config.global_tmp_dir} " +
    "INPUT=#{sam_file} " +
    "OUTPUT=#{bam_file} " +
    "VALIDATION_STRINGENCY=#{picard_validation}"
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
    "SORT_ORDER=coordinate VALIDATION_STRINGENCY=#{picard_validation}"
  end

  def dups
    "#{@config.global_java_vm} -jar -Xmx#{java_vm_mem_dups} " +
    "#{picardjars}/MarkDuplicates.jar " +
    "TMP_DIR=#{@config.global_tmp_dir} " +
    "INPUT=#{bam_file_sorted} " +
    "OUTPUT=#{bam_file_sorted_dups} " +
    "METRICS_FILE='./metric_file.picard' " +
    "VERBOSITY=ERROR " +
    "VALIDATION_STRINGENCY=#{picard_validation}"
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
    "VALIDATION_STRINGENCY=#{picard_validation}"
  end

  # Regenerate the header so we have more useful information on it
  def gen_header
    cmd = "#{@config.global_java_vm} " + 
    cmd << "-jar #{dist_dir}/#{@config.header_regen_jar} "
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
    "#{@config.global_java_vm} -cp #{@config.capture_j_classpath} " +
    "-Xmx6000M CaptureStatsBAM4 " +
    "-o #{@config.capture_stats_dir}/#{root_name} -t " +
    "#{@config.capture_chip_design} " +
    "-i #{bam_file_sorted_dups} -w -d"
  end

  def bam_reads_validator
    read_val_core + " #{bam_file_sorted_dups} #{@config.global_run_dir}" +
    " ./output/#{root_name}.read_val_log.txt"
  end

  def email_success
    extra_fn  = "email_info.txt"
    email_to  = @config.success_email_to
    cat_files = if @config.global_input_MP == 1 or @config.global_bwaaln == 1
      "marked.stats.F3.txt marked.stats.R3.txt"
    else 
      "marked.stats.txt"
    end

    # Dump here information about the run you want to send by email
    File.open(extra_fn, "w") do |f|
      f.puts "\nSEA dir: #{Dir.pwd}"
      if @config.global_input_CAP == 1
        f.puts "Cap stats: " + Dir.pwd.chomp + "/" + @config.capture_stats_dir 
      end
    end
  
    cmd = "cat #{cat_files} #{extra_fn} | "
    cmd << 'mail -s \"[OK] BFast analysis completed: ' + root_name + '\"'
    cmd << " #{email_to}"
  end

  # Clean up the analysis directory 
  # rm -rf cluster_JOBS.sh go.sh reads
  def clean_up
    cleaner_script = File.dirname(__FILE__) + "/../helpers/clean_sea_dir.sh"
  end

  private

  def dist_dir
    "#{@config.global_dist_dir}"
  end

  def update_db_tool
    dist_dir + "/bin/update_db.rb"
  end

  def picard_validation; "#{@config.global_picard_validation.upcase}" ;end

  def stats_core
    "#{@config.global_java_vm} "   +
    "-jar #{dist_dir}/#{@config.stats_s_jar} " +
    " #{bam_file_sorted_dups} "
  end

  def read_val_core
    "#{@config.global_java_vm} " +
    "-jar #{dist_dir}/#{@config.countreads_bam_reads_val_jar} "
  end

  def java_vm_mem_to_bam; "#{@config.tobam_java_vm_mem}" ;end

  def java_vm_mem_dups; "#{@config.dups_java_vm_mem}"; end

  def java_vm_mem_sort; "#{@config.sort_java_vm_mem}"; end

  def merged_bam; "./output/#{root_name}.merged.bam"; end

  def list_bams_to_merge
    sm = "_sSs_"
    t  = "output/#{sm}/bfast.reported.file.#{root_name}.#{sm}.bam"
    (1..@n_splits).inject("") {|list, s| list << "INPUT=" + t.gsub(/#{sm}/, "split#{s}") + " " }
  end

  def split_dir
    "./output/#{@curr_split}"
  end
  
  def main_bin(sub_cmd)
    "#{@config.global_bfast_bin}/bfast #{sub_cmd} "
  end

  def bfast2_bin(sub_cmd)
    "#{@config.global_bfast2_bin}/bfast #{sub_cmd} "
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
