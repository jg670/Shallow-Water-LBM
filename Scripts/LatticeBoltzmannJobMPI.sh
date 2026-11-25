#!/bin/bash
#SBATCH --job-name=LBM_MPI_Scaling
#SBATCH --partition=cpu
#SBATCH --nodes=7
#SBATCH --ntasks-per-node=64
#SBATCH --time=01:30:00
#SBATCH --mem-per-cpu=2000mb
#SBATCH --output=lbm_mpi_scaling.%j.out
#SBATCH --error=lbm_mpi_scaling.%j.err

module purge
module load compiler/intel
module load mpi/impi

export EXE="$PWD/lbm_mpi"

PROCESS_COUNTS=(1 4 9 16 25 36 49 64 81 100 121 144 169 196 225 256 289 324 361 400)
X_DIMS=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20)

WIDTH=6000
HEIGHT=6000
INTERVAL_LENGTH=10
NUM_INTERVALS=10

echo "Starting Pure MPI Scaling Test..."
echo "----------------------------------------"

for i in "${!PROCESS_COUNTS[@]}"; do
    p="${PROCESS_COUNTS[$i]}"
    x_dim="${X_DIMS[$i]}"
    
    if [ "$p" -le "$SLURM_NTASKS" ]; then
        echo "Testing with $p processes (Grid: ${x_dim}x$((p / x_dim)))..."
        
        # Execute directly with mpirun, passing the arguments to the program
        mpirun -n $p $EXE $WIDTH $HEIGHT $x_dim $INTERVAL_LENGTH $NUM_INTERVALS
        
        echo "----------------------------------------"
    else
        echo "Skipping $p processes: Exceeds requested SLURM_NTASKS ($SLURM_NTASKS)."
    fi
done

echo "All tests finished."
