require 'sequel'
require 'pathname'

# Find the proper DB
#
dbs_dir = File.join(File.dirname(File.dirname(__FILE__)), 'dbs')
DB_URL = case (`uname`.chomp)
  when /Darwin/
    "sqlite://#{dbs_dir}/test_seas.db"
  else
    "sqlite://#{dbs_dir}/seas.db"
end

DB = Sequel.connect(DB_URL)

# This are the sequel models for our SQL database
#
class Sea < Sequel::Model(:seas)
  one_to_many :key_values
end

class KeyValue < Sequel::Model(:key_values)
  many_to_one :sea
end

# Some methods to help us deal with the Sea models
#
module M_helpers
  # Add a new SEA
  def self.add(o)
    Helpers::log "Trying to add SEA"
    unless Sea[:name => o.name]
      begin 
        se = Sea.create(:name => o.name)
        se.save
        Helpers::log "SEA added"
      rescue
        Helpers::log "Problems adding SEA. Bailing out.", 1
      end   
    else
      Helpers::log "SEA already exists. Bailing out", 1
    end
  end

  # Remove a SEA
  #
  def self.remove(o)
    Helpers::log "Trying to remove SEA"
    if Sea[:name => o.name]
      begin
        Sea.find({:name => o.name}).delete
        Helpers::log "SEA removed."
      rescue
        Helpers::log "Problems removing SEA. Bailing out.", 1
      end
    else
      Helpers::log "I couldn't find the SEA.", 1
    end
  end

  # Update a key/value 
  #
  def self.update(o)
    Helpers::log "Trying to update SEA"
    att = o.key.to_sym
    if Sea[:name => o.name]
      begin
        sea = Sea.find({:name => o.name})
        if sea.respond_to?(att)
          sea[att] = o.value
          sea.save
          Helpers::log "SEA updated (level1)."
        else
          kv = KeyValue.create(:key => o.key, :value => o.value)
          sea.add_key_value kv
          sea.save
          kv.save
          Helpers::log "SEA updated (level2)."
        end
      rescue
        Helpers::log "Problems updating SEA. Bailing out."
      end
    else
      Helpers::log "I couldn't find the SEA.", 1
    end
  end
end
