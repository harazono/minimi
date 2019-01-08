#! /usr/bin/env ruby
if ARGV.size < 1
  k = 100
else
  k = ARGV[0].to_i
end
base = ['A', 'C', 'G', 'T']

r = Random.new
ref_seq     = (1..k).map{ base[r.rand(base.length)] }.join
query_seq   = (1..k).map{ base[r.rand(base.length)] }.join

puts sprintf(">rand_ref_%dbases\n%s\n>rand_qry_%dbases\n%s\n", k, ref_seq, k, query_seq)

