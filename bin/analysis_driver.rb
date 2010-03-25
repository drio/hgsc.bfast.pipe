#!/usr/bin/env ruby
#
#== Synopsis 
# This script is the entry point to perform bfast analysis   
#
#== Examples
# analysis_driver.rb -run_name RUN -c CAP -t FR -q normal -a clean
#
#== Usage 
#  analysis_driver.rb [options]
#
#  For help use: analysis_driver.rb -h
#
#== Options
#   -h, --help          Displays help message
#   -v, --version       Display the version, then exit
#   -q, --quiet         Output as little as possible, overrides verbose
#   -V, --verbose       Verbose output
#   -r, --run_name      Run_name
#   -c, --c_design      capture_design
#   -t, --se_type       sequence eventy type
#   -q, --queue         cluster queue
#   -a, --action        action to perform 
#
require 'optparse' 
require 'rdoc/usage'
require 'ostruct'
require 'date'
require 'logger'

#
# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80 
#
class App
  VERSION = '0.0.1'
  
  attr_reader :options

  def initialize(arguments, stdin)
    @arguments       = arguments
    @stdin           = stdin
    @logger          = Logger.new(STDERR)

    # Set defaults
    @options         = OpenStruct.new
    @options.verbose = false
    @options.quiet   = false
    # TO DO - add additional defaults
  end

  # Parse options, check arguments, then process the command
  def run
    if parsed_options? && arguments_valid?
      log "Start at #{DateTime.now}\n\n"
      output_options if @options.verbose

      process_arguments
      process_command
      log "Finished at #{DateTime.now}"
    else
      output_usage
    end
  end
  
  protected

    def parsed_options?
      # Specify options
      opts = OptionParser.new 
      opts.on('-v', '--version')      { output_version ; exit 0 }
      opts.on('-h', '--help')         { output_help }
      opts.on('-V', '--verbose')      { @options.verbose  = true }  
      opts.on('-q', '--quiet')        { @options.quiet    = true }

      opts.on('-r', '--run_name r')   {|r| @options.run_name = r }
      opts.on('-c', '--c_design c')   {|c| @options.c_design = c }
      opts.on('-t', '--se_type  s')   {|s| @options.se_type  = s }
      opts.on('-q', '--queue    q')   {|q| @options.queue    = q }
      opts.on('-a', '--action   a')   {|a| @options.action   = a }
            
      opts.parse!(@arguments) rescue return false
      process_options
      true
    end

    # Performs post-parse processing on options
    def process_options
      #@options.verbose = false if @options.quiet
    end
    
    def output_options
      @options.marshal_dump.each {|name, val| puts "#{name} = #{val}" }
    end

    # True if required arguments were provided
    def arguments_valid?
      true if @options.run_name && @options.c_design
              @options.se_type  && @options.queue
              @options.action   
    end
    
    # Setup the arguments
    def process_arguments
      # TO DO - place in local vars, etc
    end
    
    def output_help
      output_version
      RDoc::usage() #exits app
    end
    
    def output_usage
      RDoc::usage('usage') # gets usage from comments above
    end
    
    def output_version
      puts "#{File.basename(__FILE__)} version #{VERSION}"
    end
    
    def process_command
      # TO DO - do whatever this app does
      puts "Show time .." 
      #process_standard_input # [Optional]
    end

    def process_standard_input
      input = @stdin.read      
      # TO DO - process input
      
      # [Optional]
      # @stdin.each do |line| 
      #  # TO DO - process each line
      #end
    end

    def log(msg)
      @logger.info msg.chomp if @options.verbose     
    end
end

# Create and run the application
app = App.new(ARGV, STDIN)
app.run
