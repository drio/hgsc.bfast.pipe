#!/usr/bin/env ruby19
# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80
#
# Tool for adding/remove/updating seas in the DB
#
require 'ostruct'
require 'optparse'
require 'optparse/time'
#
require 'lib/helpers'
require 'lib/models'
require 'open-uri'

# Display SEA status in different ways
def show(o, rs) # args and reg exp of valid things to show
  Helpers::log("User requested to show #{o.show} SEAs.")
  case o.show
    when "pending"
      Sea.each {|s| puts "#{s.id}, #{s.name}, #{s.priority}" if s.started.nil? }
    when "completed"
      Sea.each {|s| puts "#{s.id}, #{s.name}" if s.completed }
    when "computing"
      Sea.each {|s| puts "#{s.id}, #{s.name}" if s.started && !s.completed }
    when "all"
      Sea.each do |s| 
        puts "#{s.id}, #{s.name}, #{s.added}, #{s.started}, #{s.completed}"
      end
    else
      Helpers::log("I can't process this 'show' cmd")
  end
end

def load_csv_pending(g_url)
  entries = []
  i=0
  open(g_url) do |f|
    f.each_line("\n") do |row|
      i = i + 1
      next unless row =~ /^0/
      row.chomp!
      ar = row.split(",")
      #Helpers::log("warning: n cols != 8 (#{ar.size} - #{i})") if ar.size != 8
      entries << ar
    end
  end
  entries
end

# Load pending analysis from gdoc 
# NOTICE this is heavily dependant on the format of the gdoc (csv)
def load_pending(g_url)
  default_url =  "http://spreadsheets.google.com/" +
                 "pub?key=0AiinUSoGtvz7dHZpd1E2cjJGYzVpbkYwN2pYWXQyb2c&gid=2" +
                 "&output=csv"
  g_url = (g_url == "default") ? default_url : g_url
  entries_processed = 0

  # Iterate over all the lines
  load_csv_pending(g_url).each_with_index do |e,i|

    r_name, s_name, prj, c_design, type, pty, d_added, comments = e
    sea_name = "#{r_name}_#{s_name}".force_encoding('ASCII-8BIT').gsub(/ /,'')
    Helpers::log("+ loading (#{i}): #{sea_name}")

    # Insert the SEA
    M_helpers::add OpenStruct.new({:name  => sea_name })

    # Insert all the key values for that SEA
    {
      :project     => prj     , :priority => pty,      :sample => s_name,
      :chip_design => c_design, :notes    => comments, :added  => Time.now,
    }.each do |key, value|
      key.to_s.force_encoding('ASCII-8BIT')
      value.force_encoding('ASCII-8BIT').chomp if !value.nil? and 
                                                  value.instance_of?(String)
      Helpers::log("updating (#{i}) key: #{key} || value: #{value}")
      M_helpers::update OpenStruct.new({:name  => sea_name,
                                        :key   => key,
                                        :value => value})
    end
    entries_processed = i
  end

  Helpers::log("+ #{entries_processed} SEAs added.")
end

# Main
#
o = OpenStruct.new

OptionParser.new do |opts|
  o.valid_actions = "add, remove, update"
  o.show_values   = "pending, completed, computing, all"
  opts.on '-a', '--action=ACTION', "Action (#{o.valid_actions})" do |a|
    o.action = a
  end

  opts.on '-k', '--key=KEY', 'Key' do |k|
    o.key = k
  end

  opts.on '-v', '--value=VALUE', 'Value' do |v|
    o.value = v
  end

  opts.on '-n', '--name=NAME', 'SEA Name' do |n|
    o.name = n
  end

  opts.on '-g', '--g=GURL', 'Gdoc URL (gdoc URL)' do |g|
    o.gdoc_url = g
  end

  opts.on '-s', '--show=WHAT', "What SEAs to show (#{o.show_values})" do |s|
    o.show = s
  end

  opts.on_tail '-h', '--help', 'Print this help' do
    puts opts
    exit 0
  end
end.parse! ARGV

h = "For help use -h."
if o.action
  ra = /add|remove|update/
  Helpers::log("I need a SEA name. #{h}" , 1) if o.action =~ ra && !o.name
  Helpers::log("I need a key/value. #{h}", 1) if o.action == 'update' &&
                                                 (!o.key || !o.value)
  Helpers::log("Action: #{o.action}")
  case o.action
    when "add"         ; M_helpers::add o
    when "remove"      ; M_helpers::remove o
    when "update"      ; M_helpers::update o
    else
      Helpers::log("I cannot process this action. #{h}")
  end
elsif !o.action and o.show
  rs = /pending|completed|computing|all/
  Helpers::log("What do you want me to show?. #{h}" , 1) unless o.show =~ rs

  show o, rs
elsif !o.action and o.gdoc_url
  load_pending o.gdoc_url
else
  Helpers::log("#{h}" , 1)
end
