#! /bin/bash
#SBATCH --partition=carter-compute
#SBATCH --output=/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/bin/slurm_logs/%x.%A.out
#SBATCH --mem=60G
#SBATCH -n 32
#SBATCH -t 14-00:00:00

# USAGE: sbatch --job-name=Hangauer2017_BT474_RNA-seq_fastq-dump 1_sra_download.sh

source activate get_data
python 1_sra_download.py