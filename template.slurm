#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --output=1iter100.out
#SBATCH --error=error.txt
#SBATCH --partition=gpuA100x4
#SBATCH --account=mzu-delta-gpu
#SBATCH --job-name=threaded-test
#SBATCH --time=00:04:00      # hh:mm:ss for the job
#SBATCH --constraint="scratch"

# p is the min number of partitions (chares that do traversal)

#nsys profile --stats=true --trace=cuda -o nsys4pe 

./build/cuda/gpu_nbody -n 10000 -m DS -i 1
# ./build/cuda/gpu_nbody -n 1000 -m DS
# ./build/cuda/gpu_nbody -n 10000 -m DS

