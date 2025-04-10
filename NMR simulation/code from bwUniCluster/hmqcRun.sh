#!/bin/bash
#SBATCH --job-name=hmqc
#SBATCH --time=10:00:00
#SBATCH --array=1-3
#SBATCH --cpus-per-task=26
#SBATCH --ntasks-per-node=3
#SBATCH --output=taskLog/hmqc_%A_%a.out
#SBATCH --partition=single
#SBATCH --mail-user=mengjia.he@kit.edu	

cd /pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid/experiments/

# Load the MATLAB module
module load math/matlab/R2023b

# Run MATLAB with your script
matlab -nodesktop -nodisplay -nosplash -r "hmqcRun(${SLURM_ARRAY_TASK_ID},'glucose'); exit"