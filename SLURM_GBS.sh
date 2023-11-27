#!/bin/bash
#SBATCH -J GTA_Barley
#SBATCH -o logs_slurm/GTA_Barley-%j.out
#SBATCH -c 8
#SBATCH -p soyagen
#SBATCH --mail-type=ALL
#SBATCH --mail-user=
#SBATCH --time=10-00:00
#SBATCH --mem=50G

module load sabre/1.000
module load cutadapt/2.1
module load bwa/0.7.17
module load samtools/1.8
module load vcftools/0.1.16
module load java/jdk/1.8.0_102
module load beagle/5.0
module load python/2.7
module load htslib/1.8
module load platypus/0.8.1.1

ulimit -S -n 40000

./fastgbs_V3.sh parameters_V3.txt
