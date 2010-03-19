# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

%w(yaml fileutils singleton digest/md5).each { |dep| require dep }

# Reads the bfast experiment config yaml and prepares config
class Config
  def initialize(config)
    %w(input global match local post tobam
       sort dups final header stats countreads capture success).each do |r|
      set config, r
    end
  end

  # Set all the config entries as methods for this class
  def set(config, r)
    puts "Loading config option: #{r}_options"
    config["#{r}_options"].each do |key, value|
      puts "loading suboption: #{r}_#{key}"
      instance_variable_set("@#{r}_#{key}", value)
    end
  end

  # If method missing, use the method name and return the value of the instance
  # variable with that name
  def method_missing(m, *args, &block)
    eval "@#{m}"
  end
end
