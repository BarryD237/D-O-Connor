#!/bin/bash
#
#SBATCH -p highmem
#SBATCH -n 8
#SBATCH --job-name=kallisto
#SBATCH -o kallisto.out # STDOUT
#SBATCH -e kallisto.err # STDERR

reference=/data/bdigby/D_O_Connor/reference/Homo_sapiens.GRCh38.cdna.all.fa.idx

while read x y; do

        base=$(basename $x .fastq.gz | sed 's/_.*//')

        kallisto quant -i $reference -o ../Quantification/"${base}" --bias $x $y

done < file_list.txt
