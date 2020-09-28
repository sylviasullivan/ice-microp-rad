#!/bin/ksh
#SBATCH --account=bb1131
#SBATCH --job-name=remapcon
#SBATCH --partition=compute2
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --output=/work/bb1131/b380873/tropic_obs/LOG_remapcon.run.%j.o
#SBATCH --error=/work/bb1131/b380873/tropic_obs/LOG_remapcon.run.%j.o
#SBATCH --time=00:15:00

# REMAP SIMULATION OUTPUT
basedir='/scratch/b/b380873/tropic_run2_restart/output_2017080800-2017080806/'
#cdo remapcon,targetGrid_global1.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapcon_global1.nc'
#cdo remapcon,targetGrid_global0.5.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapcon_global0.5.nc'
#cdo remapcon,targetGrid_r196x92.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapcon_r196x92.nc'
#cdo remapcon,targetGrid_global0.025.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapcon_global0.025.nc'

# REMAP EXTERNAL PARAMETERS
#basedir='/work/bb1131/b380459/TROPIC/extpar/'
#basedir2='/work/bb1131/b380873/tropic_obs/'
#cdo remapdis,targetGrid_global0.025.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap.nc' $basedir2'extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap_remapcon0.025.nc'

# REMAP USING WEIGHTS
for step in 61 62 63 64 65; do
    echo $step
    cdo remapcon,targetGrid_global0.025.nc,remap_weights_con.nc $basedir'CLCONV_3D_icon_tropic_00'$step'.nc' $basedir'CLCONV_3D_icon_tropic_00'$step'_remapcon_global0.025.nc'
done
