#!/usr/bin/env ruby

#require '../lib/to_yaml'
require 'fileutils'
require '../lib/config.rb'
require '../lib/helpers.rb'

clean_all = false
run_test = false
test_type = ""

def usage()
  puts "Usage:"
  puts "ruby19 " + __FILE__ + " action=actiontype testtype=test"
  puts "    actiontype : {clean - Removes all test directories}"
  puts "                 {run   - Executes specified testtype}" 
  puts "    testtype   : {fr    - Fragment test}"
  puts "    testtype   : {mp    - Matepair test}"
  puts "    testtype   : {cap   - capture test}"
  exit 1
end

def remove_dirs()
  puts "Removing test directories..."
  
  if File::exist?("fr") && File::directory?("fr")
    FileUtils.rm_rf("fr")
  end

  if File::exist?("mp") && File::directory?("mp")
    FileUtils.rm_rf("mp")
  end

  if File::exist?("cap") && File::directory?("cap")
    FileUtils.rm_rf("cap")
  end
  puts "Done."
  exit 0
end

def prepare_test_env(test_type)
  curr_dir = Dir.pwd
  if File::exist?(test_type) && File::directory?(test_type)
    FileUtils.rm_rf(test_type)
  end
  puts "Creating Dir : " + test_type
  FileUtils.mkdir(test_type)
  FileUtils.mkdir(test_type + "/reads")
  FileUtils.mkdir(test_type + "/track_jobs")

  if test_type.eql?("mp")
    `bzip2 -d ./data/mp/*`
  else
    `bzip2 -d ./data/mp/r3.*`
     if !File::exist?('./data/fr')
       FileUtils.mkdir('./data/fr')
     end
     csfasta_path = curr_dir + "/data/mp/r3.csfasta"
     qual_path = curr_dir + "/data/mp/r3.qual"
     FileUtils.ln_s(csfasta_path, './data/fr/r3.csfasta', :force => true)
     FileUtils.ln_s(qual_path, './data/fr/r3.qual', :force => true)
  end

  FileUtils.cd(test_type, :verbose => true)
end

# Function to modify run configuration for test
def fix_for_test(bf_config)
  c_path = Dir.pwd
  bf_config.gsub!(/\/h\/hsap.36.1.hg18\/hsap_36.1_hg18.fa/, "/t/test/test.fa")
  %w{28000 8000 4000}.each {|n| bf_config.gsub!(/#{n}/, "400") }
  %w{8g 4g}.each           {|n| bf_config.gsub!(/#{n}/, "1g") }

  bf_config.gsub!(/regen_jar:.+$/,
                  "regen_jar: #{c_path}" +
                  "/java/bam.header.creation/bam.header.creation.jar")
  bf_config.gsub!(/s_jar:.+$/,
                  "s_jar: #{c_path}/java/BAMStats/BAMStats.jar")
  bf_config.gsub!(/bam_reads_val_jar:.+$/,
                  "bam_reads_val_jar: #{c_path}" +
                  "/java/raw.bam.reads.validator/raw.bam.reads.validator.jar")
  bf_config.gsub!(/reads_per_file:.+$/, "reads_per_file: 480")
  bf_config
end

# Function to generate configuration YAML
def generate_yaml(test_type)

  # Location of test data
  data_fr = "/test/data/fr"
  data_mp = "/test/data/mp"
  chip_design =  "/stornext/snfs1/next-gen/software/hgsc/capture_designs" +
                 "/HD_exome/HD2_exome_target_region.bed.seq"

  puts "Generating YAML for test type : " + test_type

  curr_dir = Dir.pwd
  Dir.chdir("../../")
  yaml_inst = Yaml_template.new
  bf_config = yaml_inst.to_s
  bf_config.gsub!(/__RN__/        , "run_small_test")
  bf_config.gsub!(/__IMP__/       , test_type.eql?("mp") ? "1" : "0")
  bf_config.gsub!(/__ICAP__/      , test_type.eql?("cap") ? "1" : "0")
  bf_config.gsub!(/__RUN_DIR__/   , test_type.eql?("mp") ? Dir.pwd + data_mp :
                  Dir.pwd + data_fr)
  bf_config.gsub!(/__READS_DIR__/ , "./reads")
  bf_config.gsub!(/__OUTPUT_DIR__/, "./output")
  bf_config.gsub!(/__CD__/        , test_type.eql?("cap") ? chip_design : "" )

  bf_config = fix_for_test(bf_config)

  # Write the config yaml file to the test/{mp, fr, cap} directory
  Dir.chdir(curr_dir)
  File.open("bf.config.yaml", "w") {|f| f.write(bf_config)}
end

def create_lsf_jobs
  `../../bin/bfast.split.reads.rb ./bf.config.yaml`
  `../../bin/bfast.lsf.submit.rb ./bf.config.yaml`
end

def schedule_lsf_jobs
  puts "Sending pipeline jobs to cluster"
  
  if !File::exists?("cluster_JOBS.sh")
    puts "ERROR : Missing cluster_Jobs.sh file"
    puts "Cannot schedule lsf jobs"
    exit 2
  end
  `./cluster_JOBS.sh`
end

if ARGV.length < 1 || ARGV.length > 2
  usage()
end

ARGV.each do |line|
  if line.eql?('action=clean')
    clean_all = true
  elsif line.eql?('action=run')
    run_test = true
  elsif line.match(/testtype=/)
    test_type = line.gsub(/testtype=/, '')
  end
end

if !clean_all && !run_test
  usage()
end

if clean_all == true
  remove_dirs()
end

if run_test == true 
  if test_type.eql?(nil) || test_type.eql?("") ||
     (!test_type.eql?("fr") && !test_type.eql?("cap") && !test_type.eql?("mp"))
    usage()
  end
end

puts "Running Test Type : " + test_type

prepare_test_env(test_type)
generate_yaml(test_type)
create_lsf_jobs()
schedule_lsf_jobs()
