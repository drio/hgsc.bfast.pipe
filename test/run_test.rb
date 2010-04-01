#!/usr/bin/env ruby

require '../lib/to_yaml'

cleanAll = false
runTest = false

def usage()
  puts "Usage:"
  puts "ruby test_ruby.rb action=actiontype testtype=test"
  exit 1
end

if ARGV.length < 1 || ARGV.length > 2
  usage()
end

ARGV.each do |line|
  if line.casecmp('action=clean')
    cleanAll = true
  end
end

puts "Clean All : " + cleanAll.to_s
puts "Run Test : " + runTest.to_s

#y = Y_Fragment.new("run_small_test", "normal", 8,
#"/stornext/snfs4/next-gen/solid/bf.references/t/test/test.fa",
#"/stornext/snfs1/next-gen/drio-scratch/bfast_related/bf.pipeline.data.test/plain.fr.very.small.data",
#"output_id_like_test_drio",  "hsap36.1")
#puts y.dump
