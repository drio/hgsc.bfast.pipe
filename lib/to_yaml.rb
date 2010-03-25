#
# Here is a Humble attempt to create yaml data on the fly.
# Notice this is based on our logic.
#
# The idea is to map objects to sections and full yaml files
#
require 'yaml'
require 'pp'

# This will be mixed-in in our sections so we can dynamically
# reset values in the instance hash that captures all the 
# section key values
module ResolveDynamicKeys
  def method_missing(name, *args, &block)
    name         = name.to_s; name.gsub!(/=/,'')
    section_name = self.class.to_s.gsub(/YS_/, '').downcase
    iv           = eval "@#{section_name}_options" 
    iv[name]     = args[0]
  end
end

# You want a Fragment Yaml? or MP Yaml or Capture Yaml or ....
class Y_Fragment
  attr_accessor :input, :global

  def initialize(dynamic, input_MP)
    @sections = {}
    # input
    @sections['input']  = YS_Input.new
    @sections['input'].dynamic  = dynamic

    # global
    @sections['global'] = YS_Global.new
    @sections['global'].input_MP = input_MP

    # ...
  end

  def dump
    tmp = ""
    @sections.each {|k,v| tmp << YAML.dump(v) }
    tmp
  end
end

# Sections here Input, Global etc...
class YS_Input
  include ResolveDynamicKeys
  attr_accessor :input_options

  # Define your key values here... 
  def initialize
    @input_options = { 
                        'dynamic'  => nil        ,
                        'run_name' => "run_test" ,
                     }
  end
end

class YS_Global
  include ResolveDynamicKeys
  attr_accessor :global_options

  def initialize
    @global_options = { 
                        'threads'  => "1" ,
                        'input_MP' => 1   ,
                      }
  end
end

# y = Y_Fragment.new("great stuff", 1)
# puts y.dump
