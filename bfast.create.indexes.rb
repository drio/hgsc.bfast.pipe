#!/usr/bin/env ruby19

$:.unshift File.join(File.dirname(__FILE__))
%w(bfast.libs).each { |dep| require dep }

illumina_75bp_layouts = %w{
14 1111111111111111111111
14 1111101110111010100101011011111
14 1011110101101001011000011010001111111
14 10111001101001100100111101010001011111
14 11111011011101111011111111
14 111111100101001000101111101110111
14 11110101110010100010101101010111111
14 111101101011011001100000101101001011101
14 1111011010001000110101100101100110100111
14 1111010010110110101110010110111011
}

solid_50bp_layouts = %w{
14 1111111111111111111111 
14 111110100111110011111111111 
14 10111111011001100011111000111111 
14 11111111100101111000001100011111011 
14 1111111110001111110011111111 
14 111111011010011000011000110011111111 
14 11111111111110011101111111 
14 1111011000011111111001111011111 
14 11110110001011010011100101111101111 
14 1111111001000110001011100110001100011111
}

def error(msg)
	$stderr.puts msg
	exit 1
end

def bfast_convert(space, fasta_ref_file)
	if space == 0
		"bfast fasta2brg -f #{fasta_ref_file} -A 0 -t"
	elsif space == 1
		"bfast fasta2brg -f #{fasta_ref_file} -A 1 -t"
	else
		error "space mode incorrect for converting fasta file"
	end	
end

def bfast_index(hs, mask, i_number, nt_cs, fasta_ref_file, scratch_dir)
  cmd = ""
  cmd << 'bfast index '
  cmd << '-f FASTA_FILE '
  cmd << '-A NT_CS '
  cmd << '-i I_NUMBER '
  cmd << '-m MASK '
  cmd << '-w HSIZE '
  cmd << '-n 8 '
  cmd << '-d 0 '
  cmd << '-T TMP '
  cmd << '-t '
  cmd.gsub!(/NT_CS/     , nt_cs.to_s)
  cmd.gsub!(/MASK/      , mask)
  cmd.gsub!(/HSIZE/     , hs)
  cmd.gsub!(/I_NUMBER/  , i_number)
  cmd.gsub!(/FASTA_FILE/, fasta_ref_file)
  cmd.gsub!(/TMP/       , scratch_dir)
end

# Main
#
fasta_file = ARGV[0]
lsf  = LSFDealer.new("drd", "test")

# Color space / solid
# Create LSF job to convert ref CS to bin format
convert_dep = lsf.add_job("bfast.convert.ref.cs",
                           bfast_convert(1, fasta_file),
                           "rusage[mem=10000]")
lsf.blank
	
# Create the LSF jobs for the SOLiD indexes
i_number = 1
solid_50bp_layouts.each_slice(2) do |ly|
  hs, mask = ly.map {|e| e.chomp }
  cmd = bfast_index(hs, mask, i_number.to_s, 1, fasta_file, "/space1/tmp/")
  lsf.add_job("bfast.index.cs.#{mask}",
               cmd,
               "rusage[mem=30000]span[hosts=1]",
               [convert_dep])
  i_number += 1
end

# Sequence space / illumina
# Create LSF job to convert ref SS to bin format
convert_dep = lsf.add_job("bfast.convert.ref.cs",
                           bfast_convert(0, fasta_file),
                           "rusage[mem=10000]")
lsf.blank

i_number = 1
illumina_75bp_layouts.each_slice(2) do |ly|
  hs, mask = ly.map {|e| e.chomp }
  cmd = bfast_index(hs, mask, i_number.to_s, 0, fasta_file, "/space1/tmp/")
  lsf.add_job("bfast.index.nt.#{mask}",
               cmd,
               "rusage[mem=30000]span[hosts=1]",
               [convert_dep])
  i_number += 1
end

# Create the lsf_script
lsf.create_file
