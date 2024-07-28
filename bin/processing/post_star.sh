#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/bin/slurm_logs/%x.%A.%a.out
#SBATCH --mem=32G
#SBATCH -n 20
#SBATCH -t 14-00:00:00
#SBATCH --array=1-5%5

# USAGE: sbatch --job-name=Hangauer2017_BT474_RNA-seq_post_star post_star.sh

# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/pipelines

# Dir
star_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/processed/star

# Find all the files the .sam files
sam_files=$(find $star_dir -name "*.sam" | sort)
sam_file=$(echo $sam_files | cut -d " " -f $SLURM_ARRAY_TASK_ID)

# Grab prefix of the sam file (base file name and everything before Aligned.out.sam)
prefix=$(basename $sam_file | sed 's/Aligned.out.sam//')

# Convert sam to bam
cmd="samtools view -bS $sam_file > $star_dir/$prefix.bam"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Sort bam
cmd="samtools sort -@ 16 $star_dir/$prefix.bam -o $star_dir/$prefix.sorted.bam"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Index bam
cmd="samtools index $star_dir/$prefix.sorted.bam"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
