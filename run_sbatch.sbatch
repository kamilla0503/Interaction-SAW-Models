#! /bin/bash
#SBATCH --job-name="XY_LI_3D"
#SBATCH --account=proj_1575
#SBATCH --time=0-1:0
#SBATCH --output="Attempt0.out"%A-%a.out
#SBATCH --error="Attempt0_err.out"%A-%a.out
#SBATCH --array="100,200,300"
idx=$SLURM_ARRAY_TASK_ID

# Executable
srun testMC_3D $idx 0.54 data/XY_MC/Run_Start/3D/

