#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/bin/slurm_logs/%x.%A.%a.out
#SBATCH --mem=32G
#SBATCH -n 20
#SBATCH -t 14-00:00:00
#SBATCH --array=1-5%5

# USAGE: sbatch --job-name=Hangauer2017_BT474_RNA-seq_star_stats star_stats.sh

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/pipelines

# Dir
star_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/processed/star
#star_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/previous

# Find all the files the .sorted.bam files
bam_files=$(find $star_dir -name "*.sorted.bam" | sort)
bam_file=$(echo $bam_files | cut -d " " -f $SLURM_ARRAY_TASK_ID)

# Grab prefix of the bam file (base file name and everything before sorted.bam)
prefix=$(basename $bam_file | sed 's/.sorted.bam//')

# samtools idxstats
cmd="samtools idxstats $bam_file > $star_dir/$prefix.idxstats"
echo -e "Running command:\n$cmd\n"
eval $cmd

# samtools flagstat
cmd="samtools flagstat $bam_file > $star_dir/$prefix.flagstat"
echo -e "Running command:\n$cmd\n"
eval $cmd

# samtools stats
cmd="samtools stats $bam_file > $star_dir/$prefix.stats"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
