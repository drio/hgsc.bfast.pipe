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

def load_csv_pending(g_url)
  entries = []
  i=0
  open(g_url) do |f|
    f.each_line("\n") do |row|
      i = i + 1
      next unless row =~ /^0/
      row.chomp!
      ar = row.split(",")
      Helpers::log("warning: n cols != 8 (#{ar.size} - #{i})") if ar.size != 8
      entries << ar
    end
  end
  entries
end

# Load pending analysis from gdoc 
# NOTICE this is heavily dependant on the format of the gdoc (csv)
def load_pending(g_url)
  g_url ||= "http://spreadsheets.google.com/" +
            "pub?key=0AiinUSoGtvz7dHZpd1E2cjJGYzVpbkYwN2pYWXQyb2c&gid=2" +
            "&output=csv"

  load_csv_pending(g_url).each_with_index do |e,i|
    r_name, s_name, prj, c_design, type, pty, d_added, comments = e
    sea_name = "#{r_name}_#{s_name}".force_encoding('ASCII-8BIT')
    Helpers::log("+ loading (#{i}): #{sea_name}")
    M_helpers::add OpenStruct.new({:name  => sea_name })
    {
      :project     => prj     , :priority => pty,      :sample => s_name,
      :chip_design => c_design, :notes    => comments, :added  => Time.now,
    }.each do |key, value|
      key.to_s.force_encoding('ASCII-8BIT')
      value.force_encoding('ASCII-8BIT') if !value.nil? and value.instance_of?(String)
      Helpers::log("updating (#{i}) key: #{key} || value: #{value}")
      M_helpers::update OpenStruct.new({:name  => sea_name,
                                        :key   => key,
                                        :value => value})
    end
  end
  Helpers::log("+ {i} added: #{sea_name}")
end

o = OpenStruct.new

OptionParser.new do |opts|
  valid_actions = "add, remove, udpate, load_pending"
  opts.on '-a', '--action=ACTION', "Action (#{valid_actions})" do |a|
    o.action = a
  end

  opts.on '-k', '--key=KEY', 'Key' do |k|
    o.key = k
  end

  opts.on '-v', '--v=VALUE', 'Value' do |v|
    o.value = v
  end

  opts.on '-n', '--n=NAME', 'SEA Name' do |n|
    o.name = n
  end

  opts.on '-g', '--n=GURL', 'Gdoc URL (optional)' do |g|
    o.gdoc_url = g
  end

  opts.on_tail '-h', '--help', 'Print this help' do
    puts opts
    exit 0
  end
end.parse! ARGV

# Add: $0 -a add (just create a new Id + name)
#
# Remove: $0 -a remove -n SEA_name
#
# Update: $0 -a update -k KEY -v VALUE
# 1. Find the SEA
#   + if   model Sea has that method -> insert
#   + else use key/value table
#
# load_pending: load SE from csv file 
# 
h = "For help use -h."
Helpers::log("And action is required.#{h}", 1) if !o.action
ra = /add|remove|update/
Helpers::log("I need a SEA name. #{h}"    , 1) if !o.name && o.action =~ ra
Helpers::log("I need a key/value. #{h}"   , 1) if o.action == 'update' &&
                                                  (!o.key || !o.value)

Helpers::log("Action: #{o.action}")
case o.action
  when "add"         ; M_helpers::add o
  when "remove"      ; M_helpers::remove o 
  when "update"      ; M_helpers::update o
  when "load_pending"; load_pending o.gdoc_url
  else
    Helpers::log("I cannot process this action. #{h}")
end
