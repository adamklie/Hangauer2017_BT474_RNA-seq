#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/bin/slurm_logs/%x.%A.%a.out
#SBATCH --mem=32G
#SBATCH -n 20
#SBATCH -t 14-00:00:00
#SBATCH --array=1-5%5

# USAGE: sbatch --job-name=Hangauer2017_BT474_RNA-seq_star_align star_align.sh


# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/pipelines

ref_dir=/cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/star
fastq_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/fastq/SRP079968
out_dir=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/processed/star

# Make output directory
mkdir -p $out_dir

# Find all the files that end in R1.fastq.gz and don't have L00 in the name
fastq_files=$(find $fastq_dir -name "*R1.fastq.gz" | grep -v "L00" | sort)
fastq_file=$(echo $fastq_files | cut -d " " -f $SLURM_ARRAY_TASK_ID)

# Gunzip the fastq file
cmd="gunzip -k $fastq_file"
echo -e "Running command:\n$cmd\n"
eval $cmd
fastq_file=$(dirname $fastq_file)/$(basename $fastq_file .gz)

# Grab prefix of the fastq file (base file name and everything before _R1.fastq.gz)
prefix=$(basename $fastq_file | sed 's/_R1.fastq//')

# Run STAR
cmd="STAR --runThreadN 16 \
--genomeDir $ref_dir \
--readFilesIn $fastq_file \
--outFileNamePrefix $out_dir/$prefix"
echo -e "Running command:\n$cmd\n"
eval $cmd

# Date
date
