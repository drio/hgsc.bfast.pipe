#!/usr/bin/env ruby

require 'fileutils'
#require '../lib/config.rb'
require '../lib/helpers.rb'

clean_all = false
run_test = false
test_type = ""
$chip_design = "/stornext/snfs1/next-gen/software/hgsc/capture_designs" +
               "/HD_exome/HD2_exome_target_region.bed.seq"

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

# Method to clean the test directories
def remove_dirs()
  puts "Removing test directories..."
 
  if File::exist?("tmp") && File::directory?("tmp")
    FileUtils.rm_rf("tmp") 
  end
  exit 0
end

# Method to generate a valid but random run name
def generate_run_name(test_type)
  time = Time.new
  spot_name = "0312_" + time.year.to_s + time.month.to_s + "00" +
               time.day.to_s + "_1_SP" 
  library_name = "" 
  
  if test_type.eql?("mp")
    library_name = "ANG_TEST_1_1pA_0100" + rand(999999).to_s + "_1"
  else
    library_name = "ANG_TEST_1_1sA_0100" + rand(999999).to_s + "_1" 
  end
  return spot_name + "_" + library_name
end

# Copy data with suitable names to directories used for testing
def copy_data(analysis_dir, result_dir, run_name, test_type)
  if test_type.eql?("mp")
    `cp ./data/mp/* #{result_dir}`
    `bzip2 -d #{result_dir}/*`
     FileUtils.mv(result_dir + "/f3.csfasta", result_dir +
                  "/" + run_name + "_F3.csfasta")
     FileUtils.mv(result_dir + "/f3.qual", result_dir +
                  "/" + run_name + "_F3_QV.qual")
  else
    `cp ./data/mp/r3* #{result_dir}`
    `bzip2 -d #{result_dir}/*`
  end
  FileUtils.mv(result_dir + "/r3.csfasta", result_dir +
               "/" + run_name + "_R3.csfasta")
  FileUtils.mv(result_dir + "/r3.qual", result_dir +
               "/" + run_name + "_R3_QV.qual")
end

# Creates expected directory structure and modifies yaml for testing
def prepare_test_env(test_type)
  curr_dir = Dir.pwd
  time = Time.new  
  run_name = generate_run_name(test_type)
  spot_name = run_name.slice(/0312_\d+_1_SP/)
  lib_name = run_name.slice(/ANG_TEST_1_1\wA_\d+_1$/)

  analysis_dir = "./tmp/snfs4/next-gen/solid/analysis/solid0312/" +
  time.strftime("%Y") + "/" + time.strftime("%m") # + "/" + run_name

  result_dir = "./tmp/snfs4/next-gen/solid/results/solid0312/" + spot_name + 
               "/" + lib_name + "/results.F1B1/primary/reads" 

  sea_dir = analysis_dir + "/" + run_name

  # Create necessary directories
  FileUtils.mkdir_p(analysis_dir)
  FileUtils.mkdir_p(result_dir)
  FileUtils.mkdir_p("./tmp/snfs1/next-gen/solid/analysis/solid0312")
  FileUtils.mkdir_p("./tmp/snfs1/next-gen/solid/results/solid0312")
  FileUtils.mkdir_p("./tmp/snfs5/next-gen/solid/analysis/solid0312")
  FileUtils.mkdir_p("./tmp/snfs5/next-gen/solid/results/solid0312")

  # copy data to these directories
  copy_data(analysis_dir, result_dir, run_name, test_type)

  # build analysis driver command
  analysis_driver_cmd = "ruby19 ../bin/analysis_driver.rb -a sea_create " +
                        " -r " + run_name

  # for capture test, add chip design parameter                       
  if test_type.eql?("cap")
    analysis_driver_cmd = analysis_driver_cmd + " -c " + $chip_design
  end 
 
  puts analysis_driver_cmd 
  # Execute the command to start analysis_driver
  output = `#{analysis_driver_cmd}`

  # Copy data from result directory to analysis directory. This is 
  # required as creation of soft links to data from the analysis directory
  # to the result directory does not work in the test mode
  `cp #{result_dir}/* #{analysis_dir}/#{run_name}/input`

  # Modify the yaml for testing purpose
  modify_yaml_for_testing(sea_dir)
  #return output
  
  return sea_dir
end

# Method to modify the .yaml for testing. Reference is changed to test
# reference, distribution directory is changed to working copy
# and memory requirements are reduced
def modify_yaml_for_testing(yaml_path)
  puts "Modifying yaml for testing"
  yaml_file = yaml_path + "/bf.config.yaml"
  c_path = Dir.pwd
  c_path.slice!(/\/test$/)

  bf_config = File.open(yaml_file).read
  bf_config.gsub!(/\/h\/hsap.36.1.hg18\/hsap_36.1_hg18.fa/, "/t/test/test.fa")
  bf_config.gsub!(%r{h/hsap.36.1.hg18/bwaaln/hsap_36.1_hg18.fa}, 
                  "/t/test/bwaaln/test.fa")
  %w{28000 8000 4000}.each {|n| bf_config.gsub!(/#{n}/, "400") }
  %w{8g 4g}.each           {|n| bf_config.gsub!(/#{n}/, "1g") }
  
  bf_config.gsub!(/dist_dir:.+$/, "dist_dir: #{c_path}")
  bf_config.gsub!(/reads_per_file:.+$/, "reads_per_file: 480")

  # Send email only to the user that is testing this
  user=`id -u -n`.chomp
  bf_config.gsub!(/email_to:.+$/, "email_to: #{user}@bcm.edu")

  File.open(yaml_file, "w") {|f| f.write(bf_config)}
end

# Method to create a starting script "go.sh". Need a separate method to work
# around the errors in creating links from result to analysis directory
# generated in analysis driver
def create_starting_script(sea_dir)
  c_path = Dir.pwd
  c_path.slice!(/test$/)
  run_analysis_path = c_path + "/helpers/run_analysis.sh" 
  FileUtils.cd(sea_dir)
  `#{run_analysis_path} normal > ./go.sh`
  FileUtils.chmod(0755, "./go.sh")
  `sh ./go.sh`
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
sea_dir = prepare_test_env(test_type)
puts "Generating starting script and running the test"
create_starting_script(sea_dir)
