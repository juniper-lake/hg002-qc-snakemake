#!/bin/bash
#SBATCH -p compute
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --cpus-per-task 4
#SBATCH -o cluster_logs/slurm-%x-%j-%N.out

# USAGE: sbatch workflow/process_sample.slurm.sh <sample_id>

# set variables
SAMPLE_SHEET=$1
shift
WORKFLOWDIR=workflow
ACCOUNT=100humans
export TMPDIR=/scratch
export SINGULARITY_TMPDIR="$TMPDIR"
export SINGULARITY_BIND="$TMPDIR"

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
                --account=$ACCOUNT \
                --cpus-per-task={threads} \
                --output=cluster_logs/slurm-%x-%j-%N {resources.extra} " \
    --default-resources partition='compute' tmpdir=system_tmpdir threads=1 "extra=''" \
    --snakefile ${WORKFLOWDIR}/Snakefile \
    "$@"
