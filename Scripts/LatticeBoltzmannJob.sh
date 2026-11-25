#!/bin/bash
#SBATCH --job-name=LBM_Scaling
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=64
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=2000mb
#SBATCH --output=lbm_scaling.%j.out
#SBATCH --error=lbm_scaling.%j.err

echo "ID/NAME:    $SLURM_JOB_ID / $SLURM_JOB_NAME"
echo "USER:       $SLURM_JOB_USER"
echo "PARTITION:  $SLURM_JOB_PARTITION"
echo "TASKS:      $SLURM_TASKS_PER_NODE tasks/node x $SLURM_JOB_NUM_NODES nodes = $SLURM_NTASKS tasks"
echo "NODES:      $SLURM_JOB_NODELIST"

# Clean environment and load Intel and MPI modules
module purge
module load compiler/intel
module load mpi/impi

export UCX_MEM_CACHE=n
export I_MPI_OFI_PROVIDER=mlx
export FOR_COARRAY_CONFIG_FILE="$PWD/lbm.conf"
export EXE="$PWD/lbm"

# Define the sequence of process counts to loop through
PROCESS_COUNTS=(49)
CO_DIMS=(7)

# Global simulation parameters 
WIDTH=6006
HEIGHT=6006
INTERVAL_LENGTH=10
NUM_INTERVALS=10

echo "Starting Lattice Boltzmann Scaling Test..."
echo "----------------------------------------"

for i in "${!PROCESS_COUNTS[@]}"; do
    p="${PROCESS_COUNTS[$i]}"
    co_dim="${CO_DIMS[$i]}"

    if [ "$p" -le "$SLURM_NTASKS" ]; then
        echo "Testing with $p processes..."
        
        # Dynamically write the coarray config file for the current core count
        echo "-genvall -genv I_MPI_FABRICS=shm:ofi -n $p $EXE $WIDTH $HEIGHT $co_dim $INTERVAL_LENGTH $NUM_INTERVALS" > $FOR_COARRAY_CONFIG_FILE
        
        # Run the executable. 
        # Standard output (including Time and MLUPS) goes directly to the .out file.
        $EXE $WIDTH $HEIGHT $co_dim $INTERVAL_LENGTH $NUM_INTERVALS
        
        echo "----------------------------------------"
    else
        echo "Skipping $p processes: Exceeds requested SLURM_NTASKS ($SLURM_NTASKS)."
    fi
done

echo "All tests finished."
