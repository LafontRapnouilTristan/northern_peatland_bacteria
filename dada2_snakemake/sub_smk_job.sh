#!/bin/bash
#SBATCH -J snakeflow
#SBATCH -p unlimitq
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH -o snakemake_output_%j.out
#SBATCH -e snakemake_error_%j.err
#SBATCH --mail-type=ALL

# Environment
module purge
module load bioinfo/Snakemake/7.20.0
module load devel/Miniconda/Miniconda3

# Variables
CONFIG=resources_genobioinfo.yml
COMMAND="sbatch --cpus-per-task={cluster.cpus} --time={cluster.time} --mem={cluster.mem} -J {cluster.jobname} -o snake_subjob_log/{cluster.jobname}.%N.%j.out -e snake_subjob_log/{cluster.jobname}.%N.%j.err -p {cluster.partition}"
CORES=100

# Workflow
mkdir -p snake_subjob_log
snakemake -s Snakefile --use-conda --cluster-config $CONFIG --cluster "$COMMAND" --jobs $CORES --keep-going --latency-wait 50