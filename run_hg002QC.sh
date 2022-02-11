#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# set umask to avoid locking each other out of directories
umask 002

# execute snakemake
snakemake --reason -n \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --local-cores 1 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend mamba \
    --latency-wait 120 \
    --cluster "sbatch \
                --partition={resources.partition} \
                --cpus-per-task={threads} \
                --output=cluster_logs/slurm-%x-%j-%N-{rule}-{wildcards} {resources.extra} " \
    --default-resources partition='compute' tmpdir=system_tmpdir threads=1 "extra=''" \
    --snakefile workflow/hg002QC.smk
