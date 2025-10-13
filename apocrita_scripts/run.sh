#!/bin/bash
#$ -cwd
#$ -j y
#$ -pe smp 30
#$ -l h_rt=1:0:0
#$ -l h_vmem=3G
#$ -N hp_csp_run
#$ -o logs/hp_csp_run/

module load miniforge

cd HP-CSP
conda activate HP-CSP
source .venv/bin/activate

# Print job information
echo "Job started at: $(date)"
echo "Running on node: $(hostname)"
echo "Job ID: $JOB_ID"
echo "Working directory: $(pwd)"
echo "Python version: $(python --version)"

# # Set OMP threads for CPU-only execution
# export OMP_NUM_THREADS=$NSLOTS

# Run the training job
echo "Starting HP-CSP training..."

# Single line command to avoid line continuation issues
uv run main.py run-hp "HPH2PHPH2PH9P3HPHPH4P2H2PHPHPH2P3H4P4HP2H2PHP2HPHPH3P2H4PH3P2HP3HPH3PH3P2HPHPH7PHPH3PHP2HPH2PHPH2PHPH5P2HPH2PH4PHP6HPHP3H3PHP3HPH5PH6PHPHPHPH4PHP4H5P2HP2HP" -L 11 --dim 3 --time 3420 --workers 15 --viz none --snap out/fold.png

echo "Job completed at: $(date)"