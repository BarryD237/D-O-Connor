#!/bin/bash
#
#BATCH -n 8
#SBATCH --job-name=Index
#SBATCH -o index.out # STDOUT
#SBATCH -e index.err # STDERR

kallisto index Homo_sapiens.GRCh38.cdna.all.fa
