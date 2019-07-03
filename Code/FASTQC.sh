#!/bin/bash
#
#BATCH -n 8
#SBATCH --job-name=FASTQC
#SBATCH -o fastqc.out # STDOUT
#SBATCH -e fastqc.err # STDERR

for f in *fastq.gz; do

  fastqc $f ;
  
done
