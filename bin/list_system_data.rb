#!/usr/bin/env ruby19
#
# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80
# 
# Dumps to csv (STDOUT) all the info about the SEAs
# 
# 
require 'ostruct'
require 'find'
require 'lib/helpers'

# this is the list of valid keys to expect in the 
# stats files
KEYS_PER_TAG = %w{
  XX_total_reads_considered
  XX_total_reads_mapped
  XX_throughput
  XX_effective_throughput 
}.freeze

# csv header
HEADER = %w{
  name bam_path 
  F3_total_reads_considered F3_total_reads_mapped F3_throughput F3_effective_throughput
  R3_total_reads_considered R3_total_reads_mapped R3_throughput R3_effective_throughput
}

# csv delimiter
DELIMITER = ","

# Confirm if d looks like a sea_dir
# /stornext/snfs1/next-gen/solid/analysis/solid0100/
# 2010/04/0100_20100401_1_SL_ANG_NHLBI_A15165_L_1_1sA_01003281129_1
def sea_dir?(d)
  (
  d =~ %r{/\w+/\w+/[\w-]+/\w+/analysis/\w+/\d+/\d+/\w+} and
  d.split("/").size == 10
  ) ? true : false
end

# Find requires a "/" in order to traverse a dir. Weird
class String
  def extra_slash
    self[-1] == "/" ? self : "#{self}/"
  end
end

# Find bams in s_dir
def find_bams(s_dir)
  bams = []
  Find.find(s_dir) do |f|
    bams << f if File.file?(f) and
                 (f =~ %r{sorted.dups.bam$} or
                  f =~ %r{sorted.dups.with.header.bam$} or
                  f =~ %r{merged.marked.bam$})
  end
  bams
end

# Find stats files
# marked.stats.txt # marked.stats.F3.txt # marked.stats.R3.txt
def find_stats_files(s_dir)
  stats_files = []
  Dir[s_dir + "/*stats*"].each do |f|
    stats_files << f if f =~ %r{marked.stats[.F3|.R3]*.txt$}x and 
                        File.file?(f) and File.size?(f)
  end

  # Either 1 stats file
  # Or 2 stats file (R3|F3)
  if (stats_files.size == 1 and
      stats_files[0].split("/")[-1] == "marked.stats.txt") or
     (stats_files.size == 2 and
      stats_files[0].split("/")[-1] == "marked.stats.R3.txt" and
      stats_files[1].split("/")[-1] == "marked.stats.F3.txt")
    stats_files
  else
    []
  end 
end

# Parse the stats and gather key values for lims
def to_hash(files)
  data = ""
  h    = {}

  # Load files 
  files.each {|fn| data << File.open(fn).read }

  # Gets the key values
  data.scan(/^([F3|R3]\w+): ([,\w]+)$/).each do |m| 
    key, value = m
    unless KEYS_PER_TAG.include?(key.gsub(/R3|F3/, "XX"))
      return "Invalid key found while processing stats"
    end
    h[key] = value
  end

  (h.size == 4 or h.size == 8) ?
  h :
  "Not the expected # of key/values: #{h.size}. Bailing out."
end

# Dump a csv line in the proper format
def dump_csv_line(sea_dir, bams, stats)
  delimiter = ", "
  csv_line = sea_dir.split("/")[-1] + DELIMITER
  csv_line << bams[0] + DELIMITER
  tags = stats.size == 4 ? [ "F3" ] : [ "F3", "R3" ]

  tags.each do |tag|
    KEYS_PER_TAG.each do |k| 
      key = k.gsub(/XX/, tag)
      if stats[key].nil?
        return "I cannot find key: #{key}, n_keys: #{stats.size}. Bye."
      end
      csv_line << stats[key].gsub(/,/,".") + DELIMITER
    end
  end
  
  csv_line
end

# Load Dirs to use to look for SEAs
def load_volumes(a)
  volumes = []
  Helpers::log("Loading volumes")
  a.each do |d|
    if File.directory?(d)
      volumes << d
      Helpers::log("We'll look for SEAs in #{d}")
    else
      Helpers::log("Can't find dir: #{d}. Skipping")
    end
  end
  volumes
end

# Find SEA dirs in vols (volumes)
def find_sea_dirs(vols)
  seas = {}
  for v in vols
    Find.find(v.extra_slash) do |d|
      seas[d.split("/")[-1]] = OpenStruct.new({:dir => d}) if sea_dir?(d)
    end
  end
  seas
end

# Per each SEA, find the bams, the stats
# and dump a csv line with the data
def process_seas(seas)
  seas.each do |s, sd| # sea, sea_data
    Helpers::log("Working on SEA: #{s}")
    sd.bams        = find_bams(sd.dir)
    sd.stats_files = find_stats_files(sd.dir)
    Helpers::log("Found BAMs: #{sd.bams.size} STAT_FILEs: #{sd.stats_files.size}")

    # Only process SEA if ...
    if sd.bams.size == 1 and sd.stats_files.size.to_s =~ /1|2/
      h_stats = to_hash(sd.stats_files)
      if h_stats.is_a?(String)
        return s
      else
        if csv_line = dump_csv_line(sd.dir, sd.bams, h_stats)
          puts csv_line
        else
          return csv_line # This is an error
        end
      end
    else
      Helpers::log("Skipping #{s}")
    end
  end
  true
end

# Main
#
if __FILE__ == $0
  volumes = load_volumes(ARGV)
  seas = find_sea_dirs(volumes)
  puts HEADER.join(DELIMITER)
  r_value = process_seas(seas)
  Helpers::log(r_value, 1) if r_value.class == String
end
