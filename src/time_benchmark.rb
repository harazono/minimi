#! /usr/bin/env ruby
require 'benchmark'
require 'bio'

if ARGV.size < 1 then
  STDERR.print "time_benchmark.rb <test_data.fa>\n"
  STDERR.print "test_data.fa must contain at least 2 fasta sequence.\n"
  STDERR.print "time_benchmark.rb regards 1st sequence as reference and 2nd sequence as query, respectively.\n"
  exit
end

class READS
  attr_accessor :entry
  def initialize(reads_raw)
    parse(reads_raw)
  end
  def parse(reads_raw)
    @entry = Array.new
    reads_raw.each do |f|
      entry.push(f.seq)
    end
  end
end

test_reads_raw = Bio::FlatFile.open(Bio::FastaFormat, ARGV[0])
test_reads     = READS.new(test_reads_raw)
ref = test_reads.entry[0]
qry = test_reads.entry[1]
result = Benchmark.realtime do
  `~/ms_thesis/minimi/src/align #{ref} #{qry}`
end
print "#{result.round(5)}\t#{ref.size} x #{qry.size} = #{ref.size * qry.size}\n"

