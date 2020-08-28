# This function returns an array of the mean pressures associated with each model
# level between timestep1 and timestep2. Option to save this array with boolean 
# savepmean.
import sys

def meanPL(timestep1,timestep2,savepmean=False):
    import numpy as np
    import xarray as xr

    p_sum = np.zeros((120,))
    for i in np.arange(timestep1,timestep2):
        print(i)
        basedir = '/scratch/b/b380873/tropic_run5/'
        windth = xr.open_dataset(basedir + 'WINDTH_3D_icon_tropic_00' + str(i) + '.nc')
        p = windth.pres
        pmean = p.mean(dim='ncells')
        p_sum = p_sum + pmean.isel(time=0).values
    p_mean = p_sum/(timestep2 - timestep1)
    if(savepmean):
       np.savetxt('PMEAN_' + str(timestep1) + '-' + str(timestep2) + '.txt',p_mean)
    return p_mean

if __name__ == "__main__":
   if len(sys.argv) == 4:
      meanPL(int(sys.argv[1]),int(sys.argv[2]),sys.argv[3])
   elif len(sys.argv) == 3:
      meanPL(int(sys.argv[1]),int(sys.argv[2]))
   else:
      print('Incorrect number of arguments specified')
