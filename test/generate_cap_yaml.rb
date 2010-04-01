#!/usr/bin/env ruby

require '../../lib/to_yaml'

y = Y_Fr_Capture.new("run_small_test", "normal", 8,
"/stornext/snfs4/next-gen/solid/bf.references/t/test/test.fa",
"/stornext/snfs1/next-gen/drio-scratch/bfast_related/bf.pipeline.data.test/plain.fr.very.small.data",
"output_id_like_test_drio",  "hsap36.1",
"/stornext/snfs1/next-gen/software/hgsc/capture_designs/HD_exome/HD2_exome_target_region.bed.seq")
puts y.dump
