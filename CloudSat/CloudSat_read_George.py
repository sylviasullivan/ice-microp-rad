# collection of helper functions to read hdf file, e.g., from Cloudsat/Calipso
#
# author: Georgios Papavasileiou, KIT, 12 July 2017

# Calculate 2-weeks cloud incidence #
import numpy as np
from pyhdf.SD import SD, SDC
from pyhdf import HDF, VS  # if running outside of spyder,we need to manually import VS for an unclear reason
#import pprint

# read Cloudsat/Calipso cloud incidence data from a single hdf file
def read_cloudsatcalipso_hdf_file(file,var):
    # var can be 'CloudFraction' etc
    
    hdf = SD(file, SDC.READ)
    print (hdf.datasets())

    # read height and actual cloud data
    dset   = hdf.select('Height')
    height = dset[:,:]

    dset   = hdf.select(var)
    data   = np.array(dset[:,:], dtype=float)
    #pprint.pprint( dset.attributes() )

    # process cloud data according to valid range and scale factor
    data_at= dset.attributes(full=1)
    data_sf=data_at["factor"][0]
    data_vmin=data_at["valid_range"][0][0]
    data_vmax=data_at["valid_range"][0][1]
    data[data<data_vmin]=np.nan
    data[data>data_vmax]=np.nan  
    data = data / data_sf
    
    # Close dataset    
    dset.endaccess() 
    
    # read geolocation data and time
    h  = HDF.HDF(file)
    vs = h.vstart()

    xid = vs.find('Latitude')
    latid = vs.attach(xid)
    latid.setfields('Latitude')
    nrecs, _, _, _, _ = latid.inquire()
    latitude = np.array(latid.read(nRec=nrecs))
    latid.detach()
        
    lonid = vs.attach(vs.find('Longitude'))
    lonid.setfields('Longitude')
    nrecs, _, _, _, _ = lonid.inquire()
    longitude = np.array(lonid.read(nRec=nrecs))
    lonid.detach()

    timeid = vs.attach(vs.find('Profile_time')) 
    timeid.setfields('Profile_time')
    nrecs, _, _, _, _ = timeid.inquire()
    time = timeid.read(nRec=nrecs)
    timeid.detach()
    
    # Close file
    hdf.end()
    
    return data, height, longitude, latitude, time