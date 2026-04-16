#!/bin/bash
#SBATCH -J NaCl_UCl3
#SBATCH -o slurm_%j.out
#SBATCH -e slurm_%j.err
#SBATCH -N 4                    # NUMBER OF PHYSICAL NODES
#SBATCH --ntasks=8              # MPI ranks Ex) --ntasks=64       | THIS INDICATES        |
#SBATCH --cpus-per-task=64      # OpenMP        --cpus-per-task=1 | MPI = 64 / OpenMP = 1 |
#SBATCH -t 72:00:00
#SBATCH -p shared
#SBATCH --no-requeue # TO PREVENT INADVERTENT RESTART

# ===========================================
# CALCULATION DIRECTORY (ADJUST ACCORDINGLY)
# ===========================================
WORKDIR=/home/toh8/OpenMC/Calculation/OPT05_DEPLETION_FP_REMOVAL/91_NaCl_UCl3/91_01_NaCl_UCl3_DP_URR_OFF/91_01_URR_OFF_TRIAL11_CECM
NUCDATA=/home/toh8/OpenMC/nuclear-data

cd $WORKDIR
echo "[INFO] CWD = $(pwd)"

# ===========================================
# ASSIGN PROPER OpenMP THREAD NUMBERS
# ===========================================
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
export OMP_PLACES=cores
export OMP_PROC_BIND=spread
echo "[INFO] OMP_NUM_THREADS = $OMP_NUM_THREADS"

# ===========================================
# EXECUTE OPENMC CALCULATION
# ===========================================
echo "[INFO] Starting OpenMC Calculation"

srun apptainer exec \
     --pwd "$WORKDIR" \
     --bind $NUCDATA:/nuclear-data \
     /home/toh8/dagmc_build/dagmc_openmc_MPI.sif \
     python3 -u OpenMC_Depletion.py

echo "[INFO] Job finished."
