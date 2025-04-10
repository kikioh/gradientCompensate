#!/bin/bash
#SBATCH --job-name=SL_1spin
#SBATCH --time=8:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --mem-per-cpu=1500mb
#SBATCH --partition=single
#SBATCH --mail-user=mengjia.he@kit.edu	
#SBATCH --output=taskLog/SL_%A.out 	

cd /pfs/work7/workspace/scratch/ws4078-hmj2/OC/experiments


# Load the MATLAB module
module load math/matlab/R2023b

# Run MATLAB directly
matlab -nodesktop -nodisplay -nosplash -r "SL_1spin_exam(); exit" 
