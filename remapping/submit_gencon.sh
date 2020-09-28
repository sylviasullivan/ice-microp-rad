#!/bin/ksh
#SBATCH --account=bb1131
#SBATCH --job-name=gencon
#SBATCH --partition=compute2
#SBATCH --nodes=1
#SBATCH --output=/work/bb1131/b380873/tropic_vis/logs/LOG_gencon.run.%j.o
#SBATCH --error=/work/bb1131/b380873/tropic_vis/logs/LOG_gencon.run.%j.o
#SBATCH --time=00:10:00

# the gendis command is to generate distance-weighted average remap weights
basedir='/work/bb1131/b380459/TROPIC/extpar/'
#basedir='/scratch/b/b380873/tropic_run2_restart/output_2017080800-2017080806/'
cdo gencon,targetGrid_1.0.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap.nc' remap_weights_con_1.nc
