#! /usr/bin/env ruby

require "/home/harazono/ms_thesis/ms_thesis_scripts/sam_parser"

class READS
  attr_accessor :entry
  def initialize(reads_raw)
       parse(reads_raw)
  end
  def parse(reads_raw)
    @entry = Hash.new
    reads_raw.each do |f|
      @entry[f.entry_id] = f.seq
    end
  end
end


if ARGV.size < 3
  print "few args\n"
  exit
end

=begin
params = ARGV.getopts("k:")
kmer_size = params.find{|k, v| k == "k"}[1].to_i
kmer_size = 2 if kmer_size == 0
=end

samfile = File.open(ARGV[0], "r")
ref_raw = Bio::FlatFile.open(Bio::FastaFormat, ARGV[1])
ref = REF.new(ref_raw)

reads_raw = Bio::FlatFile.open(Bio::FastaFormat, ARGV[2])
reads = READS.new(reads_raw)

#puts reads.entry.find{|k, v| k.eql?("5afc8e30-75e3-4e0d-b529-f480550c55f4_Basecall_1D_template")}[0]

samfile.each do |line|
  next if line[/^@/]
  sam = SAM.new(line)
  next if sam.rname == "*" || sam.seq == "*"
  cigar         = CIGAR.new(sam.cigar)
  ref_subseq    = String.new
  ref_entireseq = String.new
  ref_entireseq = ref.refseq[sam.rname]
  unless ref_entireseq
    STDERR.print "sam.rname = \"#{sam.rname}\"\n"
    puts
    exit
  end
  c_sum = 0#whole length of cigar. including
  m_sum = 0
  i_sum = 0
  d_sum = 0
  s_sum = 0
  h_sum = 0
  i = 0
  while cigar.cigarOp[i] != nil do
    cutlen = cigar.cigarOp[i].length
    c_sum = c_sum + cutlen
    case cigar.cigarOp[i].type
      when "M" then
        m_sum = m_sum + cutlen
      when "I" then
        i_sumum = i_sum + cutlen
      when "D" then
        d_sum = d_sum + cutlen
      when "S" then
        s_sum = s_sum + cutlen
      when "H" then
        h_sum = h_sum + cutlen
      else
        break
    end
      i = i + 1
  end
  ref_subseq = ref_entireseq[sam.pos.to_i - 1, sam.pos.to_i + c_sum]
  unless ref_subseq
    STDERR.puts ref_entireseq
    STDERR.puts
    break
  end

  read = reads.entry.find{|k, v| k.eql?(sam.qname)}
  if read == nil then
    STDERR.print "Could not find read\nsam.qname = \"#{sam.qname}\"\n"
    exit
  end
  print "#{read[0]}\nCIGAR length = #{c_sum - d_sum}\nread  length = #{read[1].length}\n"


end

