# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

%w(yaml fileutils singleton digest/md5).each { |dep| require dep }

# Deals with the User Interface
class UInterface
  include Singleton

  # Loads config as argument or uses DATA
  def load_config(arguments, tool)
    @tool = tool
    (arguments.size == 1) ? File.new(arguments[0]) : error("Config not found")
  end

  def error(msg)
    puts "Error: " + msg; puts "" 
    puts usage
    exit 1
  end

  private

  def usage
    template = %Q{
    bfast hgsc pipeline
    VERSION: xversionx

    Usage: xcmdx <config_file>
    }

    template.gsub!(/xversionx/, Misc::VERSION)
    template.gsub!(/xcmdx/    , File.basename(@tool))
    template.gsub!(/^\s+/     , '')
  end
end
