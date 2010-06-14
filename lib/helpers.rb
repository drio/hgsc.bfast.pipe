# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80
#
require 'find'
require 'date'
require 'fileutils'

module Helpers

  curr_dir = FileUtils.pwd

  # TO DO: This has to by dynamic.
  # 
  SNFS             = %w(1 4 5).freeze
  L1_DIR           = `id -u -n`.chomp == "p-solid" ? "/stornext" :
                      curr_dir + "/tmp"
  SEA_DIR_TEMPLATE = "#{L1_DIR}/snfsSS/next-gen/solid/analysis/solidII"
  RAW_DIR_TEMPLATE = "#{L1_DIR}/snfsSS/next-gen/solid/results/solidII"
  SNFS_NUMBER      = "4"
  RUN_A_PATH       = `id -u -n`.chomp == "p-solid" ?
                      File.dirname($0) + "/../helpers/run_analysis.sh" :
                      curr_dir.gsub!(/test$/, "helpers/run_analysis.sh")

  def self.log(msg, bye=0)
    $stderr.puts "LOG: " + msg
    exit bye.to_i if bye != 0
  end

  # Look to see if a SEA directory for that sea already exists
  def self.dir_exists?(sea)
    found = []
    SNFS.each do |s|
      i_dir = SEA_DIR_TEMPLATE.gsub(/SS/, s).gsub(/II/, sea.instrument)
      log("I can't find #{i_dir}", 1) unless File.exists?(i_dir)
      re = %r{ ^#{i_dir}/\d+/\d+/#{sea}$ }x
      log("Looking for SEA dirs in: #{i_dir}")
      Find.find(i_dir) {|p| found << p if p =~ re and File.directory?(p) }
    end

    found
  end

  def self.find_raw_data(sea)
    found = []
    SNFS.each do |s|
      path = RAW_DIR_TEMPLATE.gsub(/SS/, s).gsub(/II/, sea.instrument)
      log("I can't find #{path}", 1) unless File.exists?(path)
      log("Looking for raw data in: #{path}")
      Find.find(path) do |path|
        found << path if File.file?(path)             and
                         !File.symlink?(path)         and
                         path =~ /(.csfasta$|.qual$)/ and
                         path =~ /reads/              and
                         sea.same_name_as?(path)
      end
    end
    found
  end

  # Get a SEA dir
  def self.a_dir_for(sea)
    SEA_DIR_TEMPLATE
    .gsub(/SS/, SNFS_NUMBER)
    .gsub(/II/, sea.instrument) + "/" +
    DateTime.now.year.to_s + "/" +
    sprintf("%.2d", DateTime.now.month) + "/" + 
    sea.to_s
  end

  # Create a bf.config.yaml and write to the proper location
  def self.dump_config(sea, bf_config)
    cfg_fname = a_dir_for(sea) + "/bf.config.yaml"
    Helpers::log "Creating config: #{cfg_fname}"
    File.open(cfg_fname, "w") {|f| f.write(bf_config)}
  end

  def self.create_dir(d)
    log("Couldn't create dir: #{d}, already exists", 1) if Dir.exists?(d)
    begin
      FileUtils.mkdir_p d
    rescue
      Helpers::log("Couldn't create dir: #{d}", 1)
    end
    log("dir created: #{d}")
  end

  def self.remove_dir(d)
    "\n### remove dir\n" + 
    "rm -rf #{d}\n\n"
  end

  def self.link_raw_data(sea_dir, raw_data)
    raw_data.each do |data_file|
      link_name = File.basename(data_file)

      # If the link is existing, bail out
      if File.exist?(sea_dir + "/" + link_name)
        Helpers::log "Link is existing. Bailing out.", 1
      end

      begin
        FileUtils.ln_s(data_file, sea_dir + "/" + link_name)
      rescue
        Helpers::log "Link is existing. Bailing out", 1
      end
    end 
  end

  def self.create_starting_script(sea)
    # TO DO
  end

  def self.kill_jobs_for(sea)
    cmd = "\n### bkill jobs for this SEA\n"
    cmd << "for i in `bjobs -w | grep #{sea.to_s} | awk '{print $1}'`\n"
    cmd << "do\n"
    cmd << "  bkill $i\n"
    cmd << "done\n"
  end

  def self.create_starting_script(sea)
    sea_dir = a_dir_for(sea)
    cmd = "\n### cd into the SEA and run the analysis\n"
    cmd << "cd #{sea_dir}" + "\n"
    cmd << RUN_A_PATH + " normal" + " > ./go.sh\n"
    cmd << "chmod 755 ./go.sh\n"
    cmd << "./go.sh\n"
  end
end
