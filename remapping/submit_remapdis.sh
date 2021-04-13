#!/bin/ksh
#SBATCH --account=bb1018
#SBATCH --job-name=remapdis
#SBATCH --partition=compute2
#SBATCH --mem=10gb
#SBATCH --nodes=1
#SBATCH --output=/work/bb1018/b380873/tropic_vis/logs/LOG_remapdis.run.%j.o
#SBATCH --error=/work/bb1018/b380873/tropic_vis/logs/LOG_remapdis.run.%j.o
#SBATCH --time=00:45:00

# REMAP SIMULATION OUTPUT
#basedir='/scratch/b/b380873/0V2M1A1R/'
basedir='/scratch/b/b380873/0V2M1A1R/'

#cdo remapdis,targetGrid_global1.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapdis_global1.nc'
#cdo remapdis,targetGrid_global0.5.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapdis_global0.5.nc'
#cdo remapdis,targetGrid_r196x92.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapdis_r196x92.nc'
#cdo remapdis,targetGrid_global0.025.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'CLCONV_3D_icon_tropic_0066.nc' $basedir'CLCONV_3D_icon_tropic_0066_remapdis_global0.025.nc'

# REMAP EXTERNAL PARAMETERS
#basedir='/work/bb1131/b380459/TROPIC/extpar/'
#basedir2='/work/bb1131/b380873/tropic_run2_output/'
#cdo remapdis,targetGrid_global0.025.nc -setgrid,icon-grid_tropic_55e170e5s40n_R2500m_gridID1.txt $basedir'extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap.nc' $basedir2'extpar_icon-grid_tropic_55e170e5s40n_R2500m_bitmap_remapdis0.025.nc'

# REMAP USING WEIGHTS
# 080706-080712, 3D = 0043..0048, 2D = 0085..0096
# 080712-080718, 3D = 0049..0054, 2D = 0097..0108
# 080718-080800, 3D = 0055..0060, 2D = 0109..0120
# 080800-080806, 3D = 0061..0066, 2D = 0121..0132
# 080806-080812, 3D = 0067..0072, 2D = 0133..0144

#name='merge'
#name='WINDTH_3D'
#name='RAD_2D'
#name='CLCONV_2D'
name='WINDTH_3D_F10MIN'
#name='RAD_3D'
res=0.025
for step in $(seq 1 18); do
    echo $step
    if [ ${#step} == 1 ]; then
         prefix='000'
    elif [ ${#step} == 2 ]; then
         prefix='00'
    elif [ ${#step} == 3 ]; then
         prefix='0'
    fi
    #cdo -remap,targetGrid_$res'.nc',remap_weights_dis_$res'.nc' $basedir$name'_icon_tropic_'$prefix$step'.nc' $basedir$name'_icon_tropic_'$prefix$step'_remapdis_'$res'.nc'
    #cdo -remap,targetGrid_$res'.nc',remap_weights_dis_$res'.nc' $basedir$name'_'$prefix$step'.nc' $basedir'UVT_icon_tropic_'$prefix$step'_remapdis_'$res'.nc'
    #cdo -remap,targetGrid_$res'.nc',remap_weights_dis_$res'.nc' -selvar,omega,pres,pres_sfc $basedir$name'_icon_tropic_'$prefix$step'.nc' $basedir'OMEGA_icon_tropic_'$prefix$step'_remapdis_'$res'.nc'
    #cdo -remap,targetGrid_0.025.nc,remap_weights_dis_0.025.nc -selvar,u,v,temp,pres,pres_sfc $basedir$name'_icon_tropic_'$prefix$step'.nc' $basedir'UVT_icon_tropic_'$prefix$step'_remapdis_0.025.nc'
    cdo -remap,targetGrid_$res'.nc',remap_weights_dis_$res'.nc' $basedir$name'_icon_tropic_'$prefix$step'.nc' $basedir$name'_icon_tropic_'$prefix$step'_'$res'deg.nc'
done
