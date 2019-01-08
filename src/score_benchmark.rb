#! /usr/bin/env ruby
require 'benchmark'
require 'bio'

if ARGV.size < 3 then
  STDERR.print "score_benchmark.rb <test_data.fa> <int kmer_size> <int seq_length(<=30)> <binary_scoretable_path>\n"
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
  def size()
    entry.size
  end
end

input_fastalist = ARGV[0]
kmer_size       = ARGV[1].to_i
seq_length      = ARGV[2].to_i - 1
binary_path     = ARGV[3] if ARGV.size >= 3

test_reads_raw = Bio::FlatFile.open(Bio::FastaFormat, input_fastalist)
test_reads     = READS.new(test_reads_raw)

cnt = 0
out_str = Array.new()

while cnt < test_reads.size do
  ref = test_reads.entry[cnt][0..seq_length]
  que = test_reads.entry[cnt + 1][0..seq_length]
  out_str.push(`~/ms_thesis/minimi/src/NWalignment #{ref} #{que} #{kmer_size} #{binary_path} `.to_i)
  cnt = cnt + 2
end
puts out_str.join(",").strip

