#!/bin/env bash

# where do you want STDOUT and STDERR to go?
#$ -o ~/maw/em-maw/sge.o.out
#$ -e ~/maw/em-maw/sge.e.out

# How much memory do you need **per core**?
#$ -l h_vmem=8G

# Number of cores you need 
#$ -pe smp 20

# Which queues would you like to submit to?
#$ -q HighMemLongterm.q,HighMemShortterm.q,LowMemLongterm.q,LowMemShortterm.q

  
# load any modules you need 
module load general/R/3.2.1
module load compilers/gcc/6.2.0
  
# Run thing that uses multiple cores 
 ~/maw/em-maw/em-maw -a DNA -i ~/maw/em-maw/hs_ref_GRCh38.p7.fa -o maws.out -k 2 -K 11 -c 1 -m 150000 -f 0 -r 1
 ~/maw/em-maw/em-maw -a DNA -i ~/maw/em-maw/hs_ref_GRCh38.p7.fa -o maws.out -k 2 -K 11 -c 0 -m 150000 -f 0 -r 1
