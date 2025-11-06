#!/bin/bash
#SBATCH --job-name=cfcwas
#SBATCH --mail-type=FAIL
#SBATCH --mail-user=your.email@example.com
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --partition=short
#SBATCH --mem=16G
#SBATCH --time=0-4:00
#SBATCH --output=logs/cfcwas_%j.out
#SBATCH --error=logs/cfcwas_%j.err

# Cell-Free Cistrome-Wide Association Studies (cfCWAS) Workflow Submission Script
# 
# This script submits the cfCWAS Snakemake workflow to a SLURM cluster.
# Update the SLURM parameters above to match your cluster configuration.

# Load required modules (adjust based on your cluster)
module load python/3.8
module load snakemake/7.0
module load conda

# Activate conda environment
source activate cfcwas
# Or if using conda activate:
# conda activate cfcwas

# Create logs directory if it doesn't exist
mkdir -p logs

# Create output directories
mkdir -p analysis/{qtl,hsq,stratAS,fusion,cwas_results,processed_bam,peaks}

# Run Snakemake workflow
# Use --cluster flag to submit jobs to SLURM
snakemake \
    --snakefile cfcwas.snakefile \
    --configfile config.yaml \
    --jobs 100 \
    --cluster "sbatch --partition={cluster.partition} --account={cluster.account} --time={cluster.time} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --mail-type=FAIL --mail-user={cluster.email}" \
    --cluster-config config.yaml \
    --latency-wait 60 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds

# Alternative: Run locally (comment out cluster submission above, uncomment below)
# snakemake \
#     --snakefile cfcwas.snakefile \
#     --configfile config.yaml \
#     --cores 4 \
#     --keep-going \
#     --printshellcmds

echo "cfCWAS workflow completed. Check logs/ directory for output."

