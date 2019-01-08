#!/usr/bin/env ruby

if ARGV.size < 1 then
  exit
end

kmer = ARGV[0].to_i
exit if kmer <= 0
exit if kmer >= 6


index_size = 5 ** kmer

print "\#ifndef _SCORE_#{kmer}MER\n\#define _SCORE_#{kmer}MER\n"
print "int score_#{kmer}mer[] = { \n"
for i in 1..index_size do
  for j in 1..index_size do
    if i == j then
       print " 1"
    else
       print "-1"
    end
    print ", " if i != index_size || j != index_size
    print "\n" if j == index_size
  end
end
print "};\n"
print "\#endif\n"

