#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# USAGE: sbatch workflow/process_sample.slurm.sh <sample_id>

WORKFLOWDIR=workflow
SAMPLE_SHEET=$1
shift

# set umask to avoid locking each other out of directories
umask 002

# execute snakemake
snakemake --config sample_sheet=${SAMPLE_SHEET} \
    --reason \
    --rerun-incomplete \
    --keep-going \
    --printshellcmds \
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
    --default-resources partition='compute' "tmpdir='/scratch'" threads=1 "extra=''" \
    --snakefile ${WORKFLOWDIR}/Snakefile \
    "$@"
