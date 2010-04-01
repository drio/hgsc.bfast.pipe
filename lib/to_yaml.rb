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

# Base class representing different run types. 
class Y_Base
  attr_accessor :input, :global

  def initialize(run_name, lsf_queue, threads, 
                 fasta_file, rn_dir, output_id, sq_type)
    @sections = {}
    # input
    @sections['input']          = YS_Input.new
    @sections['input'].run_name = run_name

    # global
    @sections['global']                 = YS_Global.new
    @sections['global'].lsf_queue       = lsf_queue
    @sections['global'].threads         = threads
    @sections['global'].fasta_file_name = fasta_file
    @sections['global'].run_dir         = rn_dir
    @sections['global'].output_id       = output_id

    @sections['post'] = YS_Post.new

    # match
    @sections['match'] = YS_Match.new

    # local
    @sections['local'] = YS_Local.new

    # tobam
    @sections['tobam'] = YS_ToBAM.new

    # sort
    @sections['sort'] = YS_Sort.new
    
    # dups
    @sections['dups'] = YS_Dup.new

    # final
    @sections['final'] = YS_Final.new
 
    #header
    @sections['header']         = YS_Header.new
    @sections['header'].sq_type = sq_type

    # stats
    @sections['stats'] = YS_Stats.new

    # countreads
    @sections['countreads'] = YS_CountReads.new

    #capture for fragment
    @sections['capture'] = YS_Capture_FR_MP.new

    # success
    @sections['success'] = YS_Success.new 
  end

  def dump
    tmp = ""
    @sections.each {|k,v| tmp << YAML.dump(v) }
    tmp
  end
end

# Class representing a fragment run type
class Y_Fragment < Y_Base
  attr_accessor :input, :global

  def initialize(run_name, lsf_queue, threads,
                 fasta_file, rn_dir, output_id, sq_type)
    super(run_name, lsf_queue, threads,
          fasta_file, rn_dir, output_id, sq_type)

    @sections['global'].input_MP  = 0
    @sections['global'].input_CAP = 0
    @sections['capture']          = YS_Capture_FR_MP.new
  end
end 

# Class representing a mate-pair run type
class Y_MP < Y_Base
  attr_accessor :input, :global

  def initialize(run_name, lsf_queue, threads,
                 fasta_file, rn_dir, output_id, sq_type)
    super(run_name, lsf_queue, threads,
          fasta_file, rn_dir, output_id, sq_type)

    @sections['global'].input_MP  = 1
    @sections['global'].input_CAP = 0
    @sections['capture']          = YS_Capture_FR_MP.new
  end
end 

# Class representing a fragment run-type with capture stats calculation
class Y_Fr_Capture < Y_Base
  attr_accessor :input, :global

  def initialize(run_name, lsf_queue, threads,
                 fasta_file, rn_dir, output_id, sq_type, chip_design)
    super(run_name, lsf_queue, threads,
          fasta_file, rn_dir, output_id, sq_type)

    @sections['global'].input_MP     = 0
    @sections['global'].input_CAP    = 1
    @sections['capture']             = YS_Capture.new
    @sections['capture'].chip_design = chip_design
  end
end 

# Sections here Input, Global etc...
class YS_Input
  include ResolveDynamicKeys
  attr_accessor :input_options

  # Define your key values here... 
  def initialize
    @input_options = { 
                        'run_name' => "run_test" ,
                     }
  end
end

class YS_Global
  include ResolveDynamicKeys
  attr_accessor :global_options

  @@pathPrefix = "/stornext/snfs1/next-gen/software/"

  def initialize
    @global_options = { 
                        'input_MP'          => 0 ,
                        'input_CAP'         => 0 ,
                        'picard_validation' => "STRICT" ,
                        'lsf_queue'         => nil ,
                        'bfast_bin'         => getBfastBinary() ,
                        'samtools_bin'      => getSamtoolsBinary() , 
                        'picardjars'        => getPicardPath() ,
                        'java_vm'           => getJavaVMPath() ,
                        'trackdir'          => "./track_jobs" ,
                        'space'             => "CS" ,
                        'fasta_file_name'   => nil ,
                        'timing'            => "ON" ,
                        'logs_dir'          => "./lsf_logs" ,
                        'run_dir'           => nil ,
                        'reads_dir'         => "./reads" ,
                        'threads'           => 1 ,
                        'output_dir'        => "./output" ,
                        'tmp_dir'           => "/space1/tmp/" ,
                        'output_id'         => nil ,
                        'reads_per_file'    => 50000 ,
                        'compress_input'    => "none" ,
                        'compress_splits'   => "none" ,
                      }
  end

  def getBfastBinary()
    return "/stornext/snfs1/next-gen/drio-scratch/" +  
           "bfast_related/versions/production/bfast"
  end

  def getSamtoolsBinary()
    return @@pathPrefix + "samtools-0.1.6"
  end

  def getPicardPath()
    return @@pathPrefix + "picard-tools/current"
  end

  def getJavaVMPath()
    return @@pathPrefix + "jdk1.6.0_01/bin/java"
  end
end

class YS_Match
  attr_accessor :match_options

  def initialize
    @match_options = {
                       'threads'       => 8 ,
                       'lsf_resources' => "rusage[mem=280]" ,
                     }
  end
end

class YS_Local
  attr_accessor :local_options

  def initialize
    @local_options = {
                       'threads'       => 8 ,
                       'lsf_resources' => "rusage[mem=280]" ,
                     }
  end
end

class YS_Post
  include ResolveDynamicKeys
  attr_accessor :post_options

  def initialize
    @post_options = {
                      'algorithm'     => 4 ,
                      'lsf_resources' => "rusage[mem=400]" ,
                      'rg_id'         => "TODO_RGID" ,
                      'rg_pl'         => "SOLiD" ,
                      'rg_pu'         => "TODO_PU" ,
                      'rg_lb'         => "TODO_RG_LIB" ,
                      'rg_ds'         => "TODO_RG_DESC" ,
                      'rg_dt'         => getCurrentTime() ,
                      'rg_sm'         => "TODO_RG_SM" ,
                      'rg_cn'         => "Baylor" ,
                    }
  end
 
  def getCurrentTime()
    time = Time.new
    formatTime = time.strftime("%Y-%m-%dT%H:%M:%S") +
    (time.utc_offset / 3600).to_s
    return formatTime
  end
end

class YS_ToBAM
  attr_accessor :tobam_options

  def initialize
    @tobam_options = {
                      'lsf_resources' => "rusage[mem=400]" ,
                      'java_vm_mem'   => "1g" ,
                    }
  end
end

class YS_Sort
  attr_accessor :sort_options

  def initialize
    @sort_options = {
                      'lsf_resources' => "rusage[mem=400]" ,
                      'java_vm_mem'   => "1g" ,
                    }
  end
end

class YS_Dup
  attr_accessor :dups_options
 
  def initialize
    @dups_options = {
                      'lsf_resources' => "rusage[mem=400]" ,
                      'java_vm_mem'   => "1g"
                    }
  end
end

class YS_Final
  include ResolveDynamicKeys
  attr_accessor :final_options
 
  def initialize
    @final_options = {
                       'lsf_resources' => "rusage[mem=400]" ,
                     }
  end
end

class YS_Header
  include ResolveDynamicKeys
  attr_accessor :header_options

  def initialize
    @header_options = {
                        'regen_jar'      => buildJarName() ,
                        'lsf_resources'  => "rusage[mem=400]" ,
                      }
  end

  def buildJarName
    path = "/stornext/snfs1/next-gen/solid/hgsc.solid.pipeline/bam.header.creation/"
    jar = "bam.header.creation.jar"
    return path + jar
  end
end

class YS_Stats
  attr_accessor :stats_options
 
  def initialize
    @stats_options = {
                       'lsf_resources' => "rusage[mem=400]" ,
                       's_jar'         => buildJarName() ,
                     }
  end

  def buildJarName()
    jarPath = "/stornext/snfs1/next-gen/solid/hgsc.solid.pipeline/BAMStats/"
    jarName = "BAMStats.jar"
    return jarPath + jarName
  end
end
  
class YS_CountReads
  attr_accessor :countreads_options

  def initialize
    buildJarName()
    @countreads_options = {
                            'bam_reads_val_jar' => buildJarName() ,
                          }
  end

  def buildJarName()
    jarPath = "/stornext/snfs1/next-gen/solid/hgsc.solid.pipeline/raw.bam.reads.validator/"
    jarName = "raw.bam.reads.validator.jar"
    return jarPath + jarName
  end
end

# Actual capture class to be used for capture tests.
class YS_Capture
  include ResolveDynamicKeys
  attr_accessor :capture_options

  def initialize
  @capture_options = {
                       'j_classpath'    => buildClassPath() ,
                       'chip_design'    => nil ,
                       'lsf_resources'  => "rusage[mem=400]" ,
                       'stats_dir'      => "cap_stats" ,
                     }
  end

  def buildClassPath()
    samPath = "/stornext/snfs1/next-gen/software/hgsc/capture_stats/sam-1.07.jar"
    picPath = "/stornext/snfs1/next-gen/software/hgsc/capture_stats/picard-1.07.jar"
    capPath = "/stornext/snfs1/next-gen/software/hgsc/capture_stats"
    return samPath + ":" + picPath + ":"  + capPath + ":."
 end
end

# Class to be used for capture options when capture is not actually required
# under fragment or mate pair scenarios. This class required to prevent
# config.rb from complaining that it did not find capture_options.
class YS_Capture_FR_MP
  attr_accessor :capture_options
 
  def initialize
    @capture_options = {
                         'lsf_resources' => "rusage[mem=400]" ,
                       }
  end
end

# Class representing success option
# We specify the email addresses to send success emails to
class YS_Success
  attr_accessor :success_options

  def initialize
  @success_options = {
                       'email_to' => "deiros@bcm.edu niravs@bcm.edu pellon@bcm.edu" ,
                     }
  end
end
