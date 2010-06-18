#!/usr/bin/env ruby19
#
# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80 
# 
# Main entry point to perform bfast analysis 
#
require 'optparse' 
require 'ostruct'
require 'date'
require 'logger'

$: << File.join(File.dirname(File.dirname($0)), "lib")
require 'load_libs'

#
class App
  VERSION = '0.0.1'
  
  attr_reader :options

  def initialize(arguments, stdin)
    @arguments     = arguments
    @stdin         = stdin
    @valid_actions = /(create|remove)/

    # Set defaults
    @options         = OpenStruct.new
    # TO DO - add additional defaults
  end

  # Parse options, check arguments, then process the command
  def run
    if parsed_options? && arguments_valid?
      log "Start at #{DateTime.now}\n"
      output_options

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

      opts.on('-r', '--run_name r')   {|r| @options.run_name = r }
      opts.on('-c', '--c_design c')   {|c| @options.c_design = c }
      opts.on('-q', '--queue    q')   {|q| @options.queue    = q }
      opts.on('-a', '--action   a')   {|a| @options.action   = a }
      opts.on('-f', '--force_mp')     { @options.force_mp = true }
            
      log "Processing arguments"
      opts.parse!(@arguments) rescue return false
      log "Parsing options"
      process_options
      true
    end

    # Performs post-parse processing on options
    def process_options
    end
    
    def output_options
      @options.marshal_dump.each {|name, val| log "#{name} = #{val}" }
    end

    # True if required arguments were provided
    def arguments_valid?
      true
    end

    # Place arguments in instance variables
    def process_arguments
      @r_name   = @options.run_name
      @sea      = @r_name.nil? ? nil : Sequence_event.new(@r_name)
      @c_design = @options.c_design || nil
      @queue    = @options.queue    || "normal"
      @action   = @options.action
      @force_pe = @options.force_mp || false
      log "Forcing MP mode detected" if @force_pe
    end
    
    def output_help
      output_version
      RDoc::usage() #exits app
    end
    
    def output_usage
      puts DATA.read
    end
    
    def output_version
      puts "#{File.basename(__FILE__)} version #{VERSION}"
    end
    
    def process_command
      error "Not valid action" unless @action =~ @valid_actions
      Driver_actions.new(@action).get_action.run(params_to_hash)
    end

    def process_standard_input
      input = @stdin.read      
      # TO DO - process input
      
      # @stdin.each do |line| 
      #  # TO DO - process each line
      #end
    end

    def params_to_hash
      {
        :r_name   => @r_name  ,
        :c_design => @c_design,
        :queue    => @queue   ,
        :action   => @action  ,
        :sea      => @sea     ,
        :force_mp => @force_pe,
      }   
    end

    def log(msg)
      Helpers::log msg.chomp
    end

    def error(msg)
      $stderr.puts "ERROR: " + msg + "\n\n"; output_usage; exit 1
    end
end

# Create and run the application
app = App.new(ARGV, STDIN)
app.run

__END__
Usage: 
  analysis_driver.rb [options]

Options:
 -h, --help          Displays help message
 -v, --version       Display the version, then exit

 -r, --run_name      Run_name
 -a, --action        action to perform 

 -f, --force_pe      Force MP despite the SE is a PE
 -c, --c_design      capture_design [optional]
 -q, --queue         cluster queue  [def: normal]

Valid actions:
 create: create the analysis dir and config file
         $ analysis_driver.rb -a sea_create -r RUN 
         $ analysis_driver.rb -a sea_create -r RUN -c C_DESIGN_DIR

 remove: check if analysis exists
         $ analysis_driver.rb -a sea_remove -r RUN
