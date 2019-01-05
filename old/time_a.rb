#! ruby -Ku -w
require 'benchmark'
require 'fileutils'

def format_duration s
  si = s.to_i
  h = si / 3600
  m = si % 3600 / 60
  (h > 0 ? "#{h}hour" : "") + (m > 0 ? "#{m}min" : "") + "#{si % 60}sec"
end

def bench
  Benchmark.bm 10 do |r|
    r.report "cfd" do
      system "a.exe"
    end
  end
end

def reduce_data
  Dir["*.plt"]
end

def main
  dir = File.readlines("in2d.txt")[30][/\d+\s+\K\w+/]
  dir &&= dir.strip
  if dir
    puts "dir: #{dir}"
    FileUtils.mkdir_p dir
    files = ["in2d.txt", "grid.csv", "grid.png"]
    files.each { |f| FileUtils.cp f, dir if File.file?(f) }
  end
  start  = Time.now
  system "a.exe"
  finish = Time.now
  elap = finish - start

  puts "start : #{start}"
  puts "finish: #{finish}"
  puts "elap: %ds (%s)" % [elap, format_duration(elap)]
end

main

__END__
