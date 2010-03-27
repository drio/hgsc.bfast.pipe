require 'singleton'

class Logger
  include Singleton

  def initialize
    @log = $stdin
  end

  def log(msg)
    @log.puts msg
  end
end

