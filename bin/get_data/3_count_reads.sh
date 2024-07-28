# mamba activate get_data

seqkit stats \
--threads 24 \
--tabular /cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/fastq/SRP079968/*/*.fastq.gz \
> /cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq/fastq/SRP079968/fastq_stats.tsv