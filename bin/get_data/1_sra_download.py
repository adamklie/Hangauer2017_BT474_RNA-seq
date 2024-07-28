import os
import sys
import glob
import subprocess
from tqdm.autonotebook import tqdm
from pysradb.sraweb import SRAweb

# Choose the current dataset we are working with
dataset_name = "Hangauer2017_BT474_RNA-seq"
srp_id = "SRP079968"

# Set-up directories
base_dir = "/cellar/users/aklie/data/datasets/Hangauer2017_BT474_RNA-seq"
cwd = os.path.join(base_dir, "bin/get_data")
fastq_dir = os.path.join(base_dir, "fastq")

# Connect to SRA
db = SRAweb()

# Grab the metadata for the SRP
metadata = db.sra_metadata(srp_id, detailed=True)

# # Download non-diabetic samples (`sra` files) to start
db.download(df=metadata, out_dir=fastq_dir, skip_confirmation=True)

# Parameters for download
tmp_dir = "/cellar/users/aklie/tmp/fastq-dump"
gzip = True
split_files = True
threads = 32

# Loop through and fastqdump out each SRA download file within the subdirectories of the fastq_dir
for sra_file in glob.glob(os.path.join(fastq_dir, srp_id, "*", "*.sra")):
    sra_dir = os.path.dirname(sra_file)
    if gzip:
        cmd = f"parallel-fastq-dump --threads {threads} --outdir {sra_dir} --split-files --tmpdir {tmp_dir} --gzip -s {sra_file}"
    else:
        cmd = f"parallel-fastq-dump --threads {threads} --outdir {sra_dir} --split-files --tmpdir {tmp_dir} -s {sra_file}"
    print(cmd)
    subprocess.run(cmd, shell=True)
