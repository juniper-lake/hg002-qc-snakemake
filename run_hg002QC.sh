#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# set umask to avoid locking each other out of directories
umask 002

LOCKFILE=conditions/hg002QC.lock

# add lockfile to directory to prevent multiple simultaneous jobs
lockfile -r 0 ${LOCKFILE} || exit 1
trap "rm -f ${LOCKFILE}; exit" SIGINT SIGTERM ERR EXIT

# execute snakemake
snakemake --reason \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
    --nolock \
    --local-cores 4 \
    --jobs 500 \
    --max-jobs-per-second 1 \
    --use-conda --conda-frontend mamba \
    --use-singularity --singularity-args '--nv '\
    --latency-wait 120 \
    --cluster "sbatch \
                --partition={resources.partition} \
                --cpus-per-task={threads} \
                --output=cluster_logs/slurm-%x-%j-%N {resources.extra} " \
    --default-resources partition='compute' tmpdir=system_tmpdir threads=1 "extra=''" \
    --snakefile workflow/hg002QC.smk "$@"
