#!/bin/bash
#SBATCH --job-name=hsqc
#SBATCH --time=3:00:00
#SBATCH --hint=nomultithread
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=1
#SBATCH --ntasks=3
#SBATCH --mem-per-cpu=1950mb
#SBATCH --partition=multiple_il
#SBATCH --mail-user=mengjia.he@kit.edu	
#SBATCH --output=taskLog/hsqc_%A.out 

cd /pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid/experiments/

OUTPUT_PATH="/pfs/work7/workspace/scratch/ws4078-hmj2/NMR_liquid/taskLog"

# Load the MATLAB module
module load math/matlab/R2023b


for SLURM_PROCID in $(seq 0 $(($SLURM_NTASKS - 1))); do

    
    TASK_ID=$((SLURM_PROCID + 1))
    OUTPUT_FILE="${OUTPUT_PATH}/hsqc_${SLURM_JOB_ID}_${TASK_ID}.out"
    
    
    
    # Run MATLAB with your script
    srun --ntasks=1 --nodes=1 --cpus-per-task=64\
         matlab -nodesktop -nodisplay -nosplash -r "hsqcRun($TASK_ID,'glucose'); exit" > $OUTPUT_FILE &
         
    echo "Running task $TASK_ID on node $SLURMD_NODENAME"
done

wait