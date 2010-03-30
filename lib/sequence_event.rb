#!/usr/bin/ruby

# This class encapsulates a Sequence Event name
#
# Example: 0312_20100211_1_SP_ANG_LVNC109718_1_1sA_01003280944_3
#
# constructor: string with the name of the sequence event
#
# public methods:
#
#  + instrument     : return instrument number (0312)
#  + l_barcode      : returns lims barcode (01003280944_3)
#  + spot?          : returns true if seq_event is a spot, false otherwise 
#  + slide?         : returns true if seq_event is slide, false otherwise 
#  + fr?            : returns true if seq_event is fragment, false otherwise 
#  + mp?            : returns true if seq_event is mate pair, false otherwise 
#  + to_s           : returns the run name for this sequence event

class Sequence_event
  def initialize(seq_event)
    @run_name = seq_event
  end

  def instrument
    return @run_name.slice(/^\d+/)
  end

  def l_barcode
    return @run_name.slice(/\d+_\d+$/)
  end

  def spot?
    if @run_name.match(/^\d+_\d+_\d+_SP_/)
      return true
    else
      return false
    end
  end

  def slide?
    if @run_name.match(/^\d+_\d+_\d+_SL_/)
      return true
    else
      return false
    end
  end

  def fr?
    if @run_name.match(/sA_\d+_\d$/)
      return true
    else
      return false
    end
  end

  def mp?
    if @run_name.match(/pA_\d+_\d$/) 
      return true
    else
      return false
    end
  end

  def to_s
    return @run_name
  end
end
