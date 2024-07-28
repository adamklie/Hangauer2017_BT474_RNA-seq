#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/bin/slurm_logs/%x.%A.out
#SBATCH --mem=32G
#SBATCH -n 20
#SBATCH -t 14-00:00:00

# USAGE: sbatch --job-name=Hangauer2017_BT474_RNA-seq_generate_star_index generate_star_index.sh


# Date
date
echo -e "Job ID: $SLURM_JOB_ID\n"

# Configuring env (choose either singularity or conda)
source activate /cellar/users/aklie/opt/miniconda3/envs/pipelines

cmd="STAR --runThreadN 16 --runMode genomeGenerate \
--genomeDir /cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/star \
--genomeFastaFiles /cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/allchrom.fa \
--sjdbGTFfile /cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/gencode.v19.annotation.gtf \
--sjdbOverhang 49 \
--outFileNamePrefix /cellar/users/aklie/data/ref/genomes/hg19/2024-mstp-bootcamp/star"
echo $cmd
eval $cmd

# Date
date
