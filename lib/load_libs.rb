# vim: set filetype=ruby expandtab tabstop=2 shiftwidth=2 tw=80

dir = File.join(File.dirname(File.dirname($0)), "lib")
puts "Lib dir: #{dir}"
Dir[File.join(dir, "*.rb")].each do |file| 
  puts "Loading #{file}"
  require File.basename(file)
end
