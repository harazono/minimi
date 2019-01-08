#! /bin/bash
#$ -S /bin/bash
#$ -cwd
./time_benchmark.rb r1k.fa
./time_benchmark.rb r10k.fa
./time_benchmark.rb r100k.fa
./time_benchmark.rb r1m.fa

