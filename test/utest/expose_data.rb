#!/usr/bin/env ruby19
#
# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80
#
require 'test/unit'
require 'fileutils'

main_dir = File.dirname(File.dirname(File.dirname(__FILE__)))
require main_dir + "/bin/list_system_data.rb"

# /tmp/snfs1/next-gen/solid/analysis/solid0100/2010/04/
# 0100_20100401_1_SL_ANG_NHLBI_A15165_L_1_1sA_01003281129_1/marked.stats.txt
# /tmp/snfs1/next-gen/solid/analysis/solid0100/2010/04/
# 0100_20100401_1_SL_ANG_NHLBI_A15165_L_1_1sA_01003281129_1/something.bam
class MockVolumes
  attr_reader :s1_dir, :s4_dir, :run_name1, :run_name2

  def initialize
    @l1_dir = "/tmp"
    @content_stats = <<EOF.gsub(/^\s+/, '')
    Uniqueness Percentage      : 75.15%
    Effective Throughput       : 2,350,834,400 bp

    BEGIN 4 LIMS
    F3_total_reads_considered: 84,715,031
    F3_total_reads_mapped: 62,563,274
    F3_throughput: 3,128,163,700
    F3_effective_throughput: 1,000,000,000
    END 4 LIMS

    Logging Quality Distrubution of Mapped Reads to : .
    Logging Quality Distrubution of Mapped Reads With M
EOF
  end

  def clean_dirs
    FileUtils.rm_rf @l1_dir + "/snfs1"
    FileUtils.rm_rf @l1_dir + "/snfs4"
  end

  def create_snfs1
    @run_name1 = "0100_20100401_1_SL_ANG_NHLBI_A15165_L_1_1sA_01003281129_1"
    @s1_dir    = @l1_dir +
                 "/snfs1/next-gen/solid/analysis/solid0100/" +
                 "2010/04/#{@run_name1}"
    FileUtils.mkdir_p @s1_dir
  end

  def create_snfs4
    @run_name2 = "0312_20100211_1_SP_ANG_LVNC109718_1_1sA_01003280944_3"
    @s4_dir  = @l1_dir +
               "/snfs4/next-gen/solid/analysis/solid0312/2010/04/" +
               "#{@run_name2}"
    FileUtils.mkdir_p @s4_dir
  end

  def create_some_bams(vn)
    d = vol_to_use(vn) 
    FileUtils.touch d + "/some_bam.bam"
    FileUtils.mkdir_p   d + "/some_dir"
    FileUtils.touch     d + "/some_dir/another_bam.bam"
  end

  def create_valid_bams(vn)
    d = vol_to_use(vn) 
    FileUtils.touch d + "/xxxx.sorted.dups.bam"
    FileUtils.touch d + "/xxxx.sorted.dups.with.header.bam"
    FileUtils.touch d + "/merged.marked.bam"
  end

  def create_one_valid_bam(vn)
    d = vol_to_use(vn) 
    FileUtils.touch d + "/merged.marked.bam"
    d + "/merged.marked.bam"
  end

  def create_empty_stats_files(vn)
    d = vol_to_use(vn) 
    FileUtils.touch d + "/marked.stats.txt"
    FileUtils.touch d + "/marked.stats.F3.txt"
    FileUtils.touch d + "/marked.stats.R3.txt"
    FileUtils.touch d + "/marked.stats.X.txt"
  end

  def add_stats_data(vn, t="") # volume number, tag
    dump_to_stats(vol_to_use(vn), t, @content_stats)
  end

  def add_incorrect_stats(vn, t="")
    data = @content_stats.gsub(/F3_throughput/, "F3_XXXXX")
    dump_to_stats(vol_to_use(vn), t, data)
  end

  def add_only_3_stats(vn, t="")
    data = @content_stats.gsub(/F3_throughput.*$/, "")
    dump_to_stats(vol_to_use(vn), t, data)
  end

  private

  def dump_to_stats(d, t, data)
    tag = t == "" ? t : ".#{t}"
    case t
      when ""
        File.open(d + "/marked.stats#{tag}.txt", "w") {|f| f.puts data }
      when "F3"
        File.open(d + "/marked.stats#{tag}.txt", "w") {|f| f.puts data }
      when "R3"
        File.open(d + "/marked.stats#{tag}.txt", "w") do |f| 
          f.puts data.gsub(/^F3/, "R3")
        end
      else
        raise "Incorrect tag: -#{tag}- ..."
    end
  end

  def vol_to_use(vn)
    case vn
      when 1
        @s1_dir
      when 4
        @s4_dir
      else
        raise "I can't work on that volume"
        exit 1
    end
  end
end

class TestLoadVolumes < Test::Unit::TestCase
  def setup
    @m = MockVolumes.new
  end
  
  def teardown
    @m.clean_dirs
  end

  def test_no_volumes
    assert_equal([]      , load_volumes([]))
    assert_equal([]      , load_volumes(["/dont_exists"]))
    assert_equal(["/tmp"], load_volumes(["/tmp"]))
  end

  def test_find_sea_dirs
    @m.create_snfs1
    s1 = "/tmp/snfs1"
    assert_equal(1   , find_sea_dirs([s1]).size)
    assert_equal(true, find_sea_dirs([s1]).has_key?(@m.run_name1))
    assert_equal(@m.s1_dir, find_sea_dirs([s1])[@m.run_name1].dir)

    s4 = "/tmp/snfs4"
    @m.create_snfs4
    assert_equal(1, find_sea_dirs([s4]).size)
    assert_equal(2, find_sea_dirs(["/tmp"]).size)
    assert_equal(true, find_sea_dirs([s4]).has_key?(@m.run_name2))
    assert_equal(@m.s4_dir, find_sea_dirs([s4])[@m.run_name2].dir)
  end

  def test_find_bams
    @m.create_snfs1
    @m.create_some_bams(1)
    assert_equal([], find_bams(@m.s1_dir))
    @m.create_valid_bams(1)
    assert_equal(3, find_bams(@m.s1_dir).size)
  end

  def test_find_stats
    @m.create_snfs1
    @m.create_empty_stats_files(1)
    assert_equal(0, find_stats_files(@m.s1_dir).size)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_stats_data(1, "")
    assert_equal(1, find_stats_files(@m.s1_dir).size)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_stats_data(1, "R3"); @m.add_stats_data(1, "F3")
    assert_equal(2, find_stats_files(@m.s1_dir).size)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_stats_data(1, ""); @m.add_stats_data(1, "F3")
    assert_equal(0, find_stats_files(@m.s1_dir).size)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_stats_data(1, ""); @m.add_stats_data(1, "R3")
    assert_equal(0, find_stats_files(@m.s1_dir).size)
  end

  def test_extra_slash
    assert_equal("/tmp/", "/tmp".extra_slash)
  end

  def test_correct_sea_dir
    assert_equal(false, sea_dir?("/tmp"))
    d = "/stornext/snfs1/next-gen/solid/analysis/solid0100/" +
        "2010/04/0100_20100401_1_SL_ANG_NHLBI_A15165_L_1_1sA_01003281129_1"
    assert_equal(true, sea_dir?(d))
    d.gsub!(/analysis/, "results")
    assert_equal(false, sea_dir?(d))
  end

  def test_to_hash
    @m.create_snfs1
    @m.add_stats_data(1, "F3")
    s_files = [@m.s1_dir + "/marked.stats.F3.txt"]
    assert_equal(4, to_hash( s_files ).size)

    @m.add_stats_data(1, "R3")
    s_files = [@m.s1_dir + "/marked.stats.F3.txt", @m.s1_dir + "/marked.stats.R3.txt"]
    assert_equal(8, to_hash( s_files ).size)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_stats_data(1, "")
    s_files = [@m.s1_dir + "/marked.stats.txt"]
    assert_equal(4, to_hash( s_files ).size)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_incorrect_stats(1, "")
    assert_equal(String, to_hash( s_files ).class)

    @m.clean_dirs
    @m.create_snfs1
    @m.add_only_3_stats(1, "")
    assert_equal(String, to_hash( s_files ).class)
  end

  def test_dump_csv_line
    require 'csv'

    @m.create_snfs1
    bam = @m.create_one_valid_bam(1)
    @m.add_stats_data(1, "")
    s_files = [@m.s1_dir + "/marked.stats.txt"]
    csv = dump_csv_line(@m.s1_dir, find_bams(@m.s1_dir), to_hash(s_files)).split(",")

    assert_equal(6, csv.size) 
    assert_equal(@m.s1_dir.split("/")[-1], csv[0])
    assert_equal(bam, csv[1]) 
    assert_equal("84.715.031"   , csv[2]) 
    assert_equal("62.563.274"   , csv[3])
    assert_equal("3.128.163.700", csv[4])
    assert_equal("1.000.000.000", csv[5])

    @m.create_snfs1
    bam = @m.create_one_valid_bam(1)
    @m.add_stats_data(1, "F3")
    @m.add_stats_data(1, "R3")
    s_files = [@m.s1_dir + "/marked.stats.F3.txt", @m.s1_dir + "/marked.stats.R3.txt"]
    csv = dump_csv_line(@m.s1_dir, find_bams(@m.s1_dir), to_hash(s_files)).split(",")

    assert_equal(10, csv.size) 
    assert_equal(@m.s1_dir.split("/")[-1], csv[0])
    assert_equal(bam, csv[1]) 
    assert_equal("84.715.031"   , csv[2]) 
    assert_equal("62.563.274"   , csv[3])
    assert_equal("3.128.163.700", csv[4])
    assert_equal("1.000.000.000", csv[5])
    assert_equal("84.715.031"   , csv[6]) 
    assert_equal("62.563.274"   , csv[7])
    assert_equal("3.128.163.700", csv[8])
    assert_equal("1.000.000.000", csv[9])
  end
end
