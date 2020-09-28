#!/bin/ksh
#SBATCH --account=bb1131
#SBATCH --job-name=griddes
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --output=/work/bb1131/b380873/tropic_vis/logs/LOG_griddes.run.%j.o
#SBATCH --error=/work/bb1131/b380873/tropic_vis/logs/LOG_griddes.run.%j.o
#SBATCH --time=01:00:00

cdo griddes /work/bb1131/b380459/TROPIC/grids/icon-grid_tropic_55e115e5s40n_R2500m.nc >icon-grid_tropic_55e115e5s40n_R2500m_1.txt
