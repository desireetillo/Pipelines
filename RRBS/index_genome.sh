#!/usr/bin/bash

#SBATCH --mem=16g
#SBATCH --cpus-per-task=8
#SBATCH --time=12:00:00
#SBATCH --partition=ccr

module load bismark/0.22.1


# mkdir bismark_mm10/
# cd bismark_mm10/
# ln -s /fdb/igenomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa mm10.fa
# cd ../

bismark_genome_preparation --parallel 8 bismark_mm10
