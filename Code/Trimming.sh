#!/bin/bash
#
#SBATCH -n 8
#SBATCH --job-name=cutadapt
#SBATCH -o cutadapt.out # STDOUT
#SBATCH -e cutadapt.err # STDERR

while read x y; do

        base1=$(basename $x .fastq.gz)
        base2=$(basename $y .fastq.gz)

        cutadapt -a AGATCGGAAGAG -A AGATCGGAAGAG \
        -o ../trimmed_reads/"${base1}".fastq.gz \
        -p ../trimmed_reads/"${base2}".fastq.gz \
        -m 32 \
        $x $y

done < file_list.txt
