#!/usr/bin/env ruby

require '../../lib/config.rb'
require '../../lib/helpers.rb'

test_run_dir = "/stornext/snfs1/next-gen/drio-scratch/bfast_related/" +
               "bf.pipeline.data.test/plain.fr.very.small.data" 
chip_design  = "/stornext/snfs1/next-gen/software/hgsc/capture_designs" +
               "/HD_exome/HD2_exome_target_region.bed.seq"

def fix_for_test(bf_config)
  bf_config.gsub!(/\/h\/hsap.36.1.hg18\/hsap_36.1_hg18.fa/, "/t/test/test.fa")
  %w{28000 8000 4000}.each {|n| bf_config.gsub!(/#{n}/, "400") }
  %w{8g 4g}.each           {|n| bf_config.gsub!(/#{n}/, "1g") }
end

def usage()
  puts "Usage:"
  puts "ruby19 " + __FILE__ + " test_type"
  puts "  test_type : Type of test"
  puts "              Allowed Values : mp, fr, cap"
  exit 1
end

if ARGV.length != 1
  usage()
end

test_type = ARGV[0].downcase
if !test_type.eql?("fr") && !test_type.eql?("mp") && !test_type.eql?("cap")
  usage()
end

pwd = Dir.pwd
Dir.chdir("../../")
yaml_inst = Yaml_template.new
bf_config = yaml_inst.to_s
bf_config.gsub!(/__RN__/        , "run_small_test")
bf_config.gsub!(/__IMP__/       , test_type.eql?("mp") ? "1" : "0")
bf_config.gsub!(/__ICAP__/      , test_type.eql?("cap") ? "1" : "0")
bf_config.gsub!(/__RUN_DIR__/   , test_run_dir)
bf_config.gsub!(/__READS_DIR__/ , "./reads")
bf_config.gsub!(/__OUTPUT_DIR__/, "./output")
bf_config.gsub!(/__CD__/        , test_type.eql?("cap") ? chip_design : "" )

# Change the reference from human to test reference to speed up
# BAM generation. Also, change the LSF usage mode
bf_config = fix_for_test(bf_config)

# Write the config yaml file to the test/{mp, fr, cap} directory
Dir.chdir(pwd)
File.open("bf.config.yaml", "w") {|f| f.write(bf_config)}
exit 0
