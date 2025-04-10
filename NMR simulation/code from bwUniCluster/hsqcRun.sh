#!/bin/bash
#SBATCH --job-name=hsqc
#SBATCH --time=3:00:00
#SBATCH --array=1-3
#SBATCH --nodes=1
#SBATCH --cpus-per-task=40
#SBATCH --ntasks-per-node=1
#SBATCH --output=taskLog/hsqc_%A_%a.out
#SBATCH --partition=single
#SBATCH --mail-user=mengjia.he@kit.edu	

cd /pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid/experiments/

# Load the MATLAB module
module load math/matlab/R2023b

# Run MATLAB with your script
# specify the sample as glycine for channel 1
# specify the sample as glucose for channel 2
# hsqc, hsqc_grad, hsqc_grad_SL are kernel functions
# hsqc_exam is example stript with sample input
matlab -nodesktop -nodisplay -nosplash -r "hsqcRun(${SLURM_ARRAY_TASK_ID},'glucose'); exit"
