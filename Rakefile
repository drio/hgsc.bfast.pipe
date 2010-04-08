require 'sequel'

def create_db(db)
  db.create_table :seas do
    primary_key :id
    String      :name
    String      :project
    String      :priority
    String      :sample
    String      :chip_design
    Text        :notes
    DateTime    :added
    DateTime    :started
    DateTime    :completed
    String      :bam
    String      :raw_data
  end

  db.create_table :key_values do
    primary_key :id
    foreign_key :sea_id, :seas
    String      :key
    String      :value
  end
end

def log(msg)
  puts "+ LOG: #{msg}"
end

@test_db_file="dbs/test_seas.db"

desc "Remove test table and create a new one"
task :create_test_db do
  log "rm ..."
  `rm -f #{@test_db_file}`
  db = Sequel.sqlite(@test_db_file)
  log "Creating db #{@test_db_file} ..."
  create_db(db)
end
