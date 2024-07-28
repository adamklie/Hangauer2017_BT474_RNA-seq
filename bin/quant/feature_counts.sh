#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/bin/slurm_logs/%x.%A.%a.out
#SBATCH --mem=32G
#SBATCH -n 20
#SBATCH -t 14-00:00:00

# USAGE: sbatch --job-name=Hangauer2017_BT474_RNA-seq-feature_counts feature_counts.sh

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/pipelines

# Dir
star_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/processed/star
#star_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/previous

# Outdir
outdir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/results/quant/featureCounts
mkdir -p $outdir

# Get all sorted.bam files in the star_dir and output for debugging (one file per line)
bam_files=$(find $star_dir -name "*.sorted.bam" | sort)
echo -e "Bam files:\n$bam_files\n"

# run feature counts
cmd="featureCounts \
-a /cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/gencode.v19.annotation.gtf \
-G /cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/allchrom.fa \
-o $outdir/hangauer.results.counts \
$star_dir/*sorted.bam"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
