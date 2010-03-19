# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

%w(yaml fileutils singleton digest/md5).each { |dep| require dep }

# Miscelanea
module Misc
  VERSION  = "0.1"
  RG_FNAME = "./rg.txt"

  # Create a file with the @RG line that bfast will use to build the @RG type
  # in the header of the bam
  #rg_id: 0                                 ; Read group (BAM spec) 
  #rg_pl: SOLiD                             ; platform 
  #rg_pu: SOLiD0097_20081114_1_Hs_1011_MP_F3; Our uniq identifier RUN + LIMS ID
  #rg_lb: HS_1011_MP                        ; Library 
  #rg_ds: rl=25                             ; Description rl=25 OR rl=50_25
  #rg_dt: 2010-01-22T18:20:29-0600          ; date, ISO 8601
  #rg_sm: CMT-001                           ; sample / pool name
  #rg_cn: Baylor                            ; center
  def Misc.create_rg_file(config)
    rg_data = "@RG\t"
    %w{ID PL PU LB DS DT SM CN}.each do |t|
      t_value = eval("config.post_rg_" + t.downcase).to_s
      puts "#{t}:-#{t_value}--" 
      raise "Couldn't find #{t} tag in config file." unless t_value
      rg_data << "#{t}:" + t_value + "\t"
    end

    File.open(RG_FNAME, "w") {|f| f << rg_data}
  end

  # Figure out how to flag bfast to process input data
  def Misc.input_compress(config)
    case config.global_compress_input
      when "bzip2"
        "-j "
      when "gzip"
        "-z "
      else  
        " "
    end
  end

  # Figure out how to flag bfast to process output data
  def Misc.output_compress(config)
    case config.global_compress_splits
      when "bzip2"
        "-J "
      when "gzip"
        "-Z "
      else  
        " "
    end
  end

  # Generate the wildcard to find files based on the user preferences
  def Misc.wild(root, config)
    case config.global_compress_input
      when "bzip2"
        "*.#{root}.bz2 "
      when "gzip"
        "*.#{root}.gz "
      else  
        "*.#{root}"
    end
  end
end
