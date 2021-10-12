#!/bin/ksh
#SBATCH --account=bb1018
#SBATCH --job-name=Ni_norm_lt
#SBATCH --partition=compute2
#SBATCH --mem=12gb
#SBATCH --nodes=1
#SBATCH --output=/work/bb1018/b380873/tropic_vis/traj/LOG_Ni_norm_lt.%j.o
#SBATCH --error=/work/bb1018/b380873/tropic_vis/traj/LOG_Ni_norm_lt.%j.o
#SBATCH --time=08:00:00

python diff_norm_lifetime.py
