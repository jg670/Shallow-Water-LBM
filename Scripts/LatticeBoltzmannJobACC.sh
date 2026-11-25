#!/bin/bash
#SBATCH --job-name=LBM_ACC_Scaling
#SBATCH --partition=gpu_a100_short
#SBATCH --gres=gpu:4
#SBATCH --nodes=2
#SBATCH --ntasks-per-node=4
#SBATCH --time=00:20:00
#SBATCH --mem-per-gpu=93000mb
#SBATCH --output=6kGPURun/lbm_acc_scaling.%j.out
#SBATCH --error=6kGPURun/lbm_acc_scaling.%j.err

module purge
module load toolkit/nvidia-hpc-sdk/25.1-byo-compiler

unset OPAL_PREFIX
unset OMPI_MCA_opal_prefix

export MPI_ROOT=/software/bwhpc/common/mpi/openmpi/5.0.8-nvidia-25.1-nompi
export PATH=$MPI_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$MPI_ROOT/lib:$LD_LIBRARY_PATH
export CPATH=$MPI_ROOT/include:$CPATH
export PKG_CONFIG_PATH=$MPI_ROOT/lib/pkgconfig:$PKG_CONFIG_PATH

export EXE="$PWD/lbm_acc"

export UCX_MEMTYPE_CACHE=n
export UCX_IB_GPU_DIRECT_RDMA=no

PROCESS_COUNTS=(5 6 7 8)
X_DIMS=(5 6 7 8)

WIDTH=6000
HEIGHT=6000
INTERVAL_LENGTH=10
NUM_INTERVALS=10

echo "Starting OpenACC + MPI Scaling Test..."
echo "----------------------------------------"

for i in "${!PROCESS_COUNTS[@]}"; do
    p="${PROCESS_COUNTS[$i]}"
    x_dim="${X_DIMS[$i]}"
    
    if [ "$p" -le "$SLURM_NTASKS" ]; then
        echo "Testing with $p processes (Grid: ${x_dim}x$((p / x_dim)))..."
        
        # Execute directly with mpirun
        srun --mpi=pmix  -n $p $EXE $WIDTH $HEIGHT $x_dim $INTERVAL_LENGTH $NUM_INTERVALS
        
        echo "----------------------------------------"
    else
        echo "Skipping $p processes: Exceeds requested SLURM_NTASKS ($SLURM_NTASKS)."
    fi
done

echo "All tests finished."
