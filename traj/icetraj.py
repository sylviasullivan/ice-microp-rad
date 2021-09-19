import os, sys
import xarray as xr
import numpy as np
from datetime import datetime, timedelta

sys.path.append(os.path.abspath("/work/bb1018/b380873/tropic_vis/utilities/"))
from plotting_utilities import general_prefix, traj_prefix
from thermodynamic_functions import *

# Calculate IWC percentiles from Ni and ri percentiles
# ri [=] microns, Ni [=] cm-3, IWC [=] ppm
def calc_IWC( Ni_perc, ri_perc):
    rho_ice = 920  # Density of ice [kg m-3]
    ri_perc = ri_perc * 10**(-6) # Ice crystal size [m]
    Ni_perc = Ni_perc * 10**6    # Ice crystal number [m-3]
    iwc = 4/3 * np.pi * rho_ice * ri_perc**3 * Ni_perc # IWC [kg m-3]

    # Somewhat elaborate conversion to ppm
    mh2o = 18 # Molar mass of water [g mol-1]
    R = 8.31 # Gas constant [J K-1 mol-1]
    T_min = 180
    T_max = 245
    Tgrid = np.arange(T_min, T_max)

    p0 = 0.000289319
    p1 = 2.66906
    p2 = -247.724
    pmean = p0*Tgrid**p1 + p2
    factor = (pmean*10**2 / Tgrid) * mh2o / R
    #print( factor )
    #print( factor.shape )
    #print( iwc )
    iwc_ppm = iwc * 10**6# * factor
    return iwc_ppm


# Filter the input fields for where IWC is not 0 or nan and temperature is between 205 and 237 K.
def filter_iwc_temp_outflow( iwc, temperature, rhi ):
    k = np.argwhere( (~np.isnan(iwc)) )
    iwc = iwc[k[:,0]]
    temperature = temperature[k[:,0]]
    rhi = rhi[k[:,0]]
    # sylvia 13092021_changing IWC lower bound
    #j = np.argwhere((iwc > 10**(-10)) & (temperature <= 238) & (temperature >= 205))
    j = np.argwhere((iwc > 10**(-8)) & (temperature <= 238) & (temperature >= 205))
    return iwc[j[:,0]], temperature[j[:,0]], rhi[j[:,0]]


# Filter the input fields for where IWC is not 0 or nan and temperature is colder than 205 K.
def filter_iwc_temp_insitu( iwc, temperature, rhi ):
    k = np.argwhere( (~np.isnan(iwc)) )
    iwc = iwc[k[:,0]]
    temperature = temperature[k[:,0]]
    rhi = rhi[k[:,0]]
    # sylvia 13092021_changing IWC lower bound
    #j = np.argwhere((iwc > 10**(-10)) & (temperature <= 205))
    j = np.argwhere((iwc > 10**(-8)) & (temperature <= 205))
    return iwc[j[:,0]], temperature[j[:,0]], rhi[j[:,0]]


# Calculate ICNC in units of [L-1] and filter out instances where the ICNC value is not 0 or nan
# and temperature is between 205 and 237 K.
def filter_icnc_temp_outflow( icnc_in, rho, temperature, rhi ):
    icnc = icnc_in * rho / 1000.
    k = np.argwhere( (~np.isnan(icnc)) )
    icnc = icnc[k[:,0]]
    temperature = temperature[k[:,0]]
    rhi = rhi[k[:,0]]
    j = np.argwhere((icnc > 10**(-10)) & (temperature <= 238) & (temperature >= 205))
    return icnc[j[:,0]], temperature[j[:,0]], rhi[j[:,0]]


# Calculate ICNC in units of [L-1] and filter out instances where the ICNC value is not 0 or nan
# and temperature is less than 205 K.
def filter_icnc_temp_insitu( icnc_in, rho, temperature, rhi ):
    icnc = icnc_in * rho / 1000.
    k = np.argwhere( (~np.isnan(icnc)) )
    icnc = icnc[k[:,0]]
    temperature = temperature[k[:,0]]
    rhi = rhi[k[:,0]]
    j = np.argwhere((icnc > 10**(-10)) & (temperature <= 205))
    return icnc[j[:,0]], temperature[j[:,0]], rhi[j[:,0]]


# Read the ice water content along the CLaMS or ICON trajectories.
# Below it is assumed that the ice water content has units of kg kg-1. Factor of 10**6 to convert to ppmv
# clams is a boolean for whether we look at ICON or CLaMS trajectories
# outflow is a boolean for whether we use the warmer T range in the filter functions above.
def read_iwctraj( set_name, time_pt, clams, outflow ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    if clams == True:
        start_str = 'cirrus_tst'
        end_str = '_trim_extract_clams.nc'
        iwc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = 10**6*( traj['IWC_hom'] + traj['IWC_het'] + traj['IWC_pre'] ).values.flatten()
            iwc.extend( vals )
            vals = traj['T'].values.flatten()
            t.extend( vals )
            vals = traj['RHI'].values.flatten()
            rhi.extend( vals )
    else:
        start_str = 'traj_tst'
        end_str = '_trim_extract.nc'
        iwc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = 10**6*traj['qi'].values.flatten()
            iwc.extend( vals )
            vals = traj['t'].values.flatten()
            t.extend( vals )
            p = traj['p'].values.flatten()
            qv = traj['qv'].values.flatten()
            vals = calc_RHi( vals, p, qv )
            rhi.extend( vals )
    iwc = np.array( iwc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the IWC value is not 0 or nan
    if outflow == True:
       iwc, t, rhi = filter_iwc_temp_outflow( iwc, t, rhi )
    else:
       iwc, t, rhi = filter_iwc_temp_insitu( iwc, t, rhi )
    print( iwc.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(iwc)), np.nanmean(np.log10(iwc)), np.nanmax(np.log10(iwc)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, iwc


# Read the sedimentation tendency along the CLaMS or ICON trajectories.
# Below it is assumed that the sedimentation tendency has units of kg kg-1 s-1. Factor of 10**6 to convert to ppmv s-1
def read_qseditraj( set_name, time_pt, clams ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    if clams == True:
        start_str = 'cirrus_tst'
        end_str = '_trim_extract_clams.nc'
        qsedi = []
        qsedo = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = 10**6*traj['qsedi'].values.flatten()
            qsedi.extend( vals )
            vals = 10**7*traj['qsedo'].values.flatten()
            qsedo.extend( vals )
            vals = traj['T'].values.flatten()
            t.extend( vals )
            vals = traj['RHI'].values.flatten()
            rhi.extend( vals )
    else:
        start_str = 'traj_tst'
        end_str = '_trim_extract.nc'
        qsedi = []
        qsedo = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = 10**6*traj['qsedi'].values.flatten()
            qsedi.extend( vals )
            vals = 10**6*traj['qsedo'].values.flatten()
            qsedo.extend( vals )
            vals = traj['t'].values.flatten()
            t.extend( vals )
            p = traj['p'].values.flatten()
            qv = traj['qv'].values.flatten()
            vals = calc_RHi( vals, p, qv )
            rhi.extend( vals )
    qsedi = np.array( iwc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the IWC value is not 0 or nan
    iwc, t, rhi = filter_iwc_temp( iwc, t, rhi )
    print( qsedi.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(qsedi)), np.nanmean(np.log10(qsedi)), np.nanmax(np.log10(qsedi)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, qsedi


# Read the ice water content along the CLaMS or ICON trajectories only for times time0-timef.
# Below it is assumed that the ice water content has units of kg kg-1. Factor of 10**6 to convert to ppmv
# clams is a boolean for whether we look at ICON or CLaMS trajectories
# outflow is a boolean for whether we use the warmer T range in the filter functions above.
def read_iwctraj_sub( set_name, time_pt, clams, outflow, time0, timef ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    if clams == True:
        start_str = 'cirrus_tst'
        end_str = '_trim_extract_clams.nc'
        iwc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = ( traj['IWC_hom'] + traj['IWC_het'] + traj['IWC_pre'] ).sel( time=slice(time0, timef) ).values.flatten()
            iwc.extend( vals*10**6 )
            vals = traj['T'].sel( time=slice(time0, timef) ).values.flatten()
            t.extend( vals )
            vals = traj['RHI'].sel( time=slice(time0, timef) ).values.flatten()
            rhi.extend( vals )
    else:
        start_str = 'traj_tst'
        end_str = '_trim_extract_dt.nc'
        iwc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = traj['qi'].sel( time=slice(time0, timef) ).values.flatten()
            iwc.extend( vals*10**6 )
            vals = traj['t'].sel( time=slice(time0, timef) ).values.flatten()
            t.extend( vals )
            p = traj['p'].sel( time=slice(time0, timef) ).values.flatten()
            qv = traj['qv'].sel( time=slice(time0, timef) ).values.flatten()
            vals = calc_RHi( vals, p, qv )
            rhi.extend( vals )
    iwc = np.array( iwc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the IWC value is not 0 or nan
    if outflow == True:
       iwc, t, rhi = filter_iwc_temp_outflow( iwc, t, rhi )
    else:
       iwc, t, rhi = filter_iwc_temp_insitu( iwc, t, rhi )
    print( iwc.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(iwc)), np.nanmean(np.log10(iwc)), np.nanmax(np.log10(iwc)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, iwc


# Read the ice water content along the full set of ICON trajectories.
# Below it is assumed that the ice water content has units of kg kg-1. Factor of 10**6 to convert to ppmv
# outflow is a boolean for whether we use the warmer T range in the filter functions above.
def read_iwctraj_all( set_name, time_pt, outflow ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    start_str = 'traj_tst'
    end_str = '.nc'
    iwc = []
    t = []
    rhi = []
    for file_id in np.arange(1,6):
        filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
        traj = xr.open_dataset( basedir + filename )
        vals = 10**6*traj['qi'].values.flatten()
        iwc.extend( vals )
        vals = traj['t'].values.flatten()
        t.extend( vals )
        p = traj['p'].values.flatten()
        qv = traj['qv'].values.flatten()
        vals = calc_RHi( vals, p, qv )
        rhi.extend( vals )
    iwc = np.array( iwc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the IWC value is not 0 or nan
    if outflow == True:
       iwc, t, rhi = filter_iwc_temp_outflow( iwc, t, rhi )
    else:
       iwc, t, rhi = filter_iwc_temp_insitu( iwc, t, rhi )
    print( iwc.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(iwc)), np.nanmean(np.log10(iwc)), np.nanmax(np.log10(iwc)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, iwc


# Read the ice crystal number concentrations along the CLaMS or ICON trajectories
def read_icnctraj( set_name, time_pt, clams ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    if clams == True:
        start_str = 'cirrus_tst'
        end_str = '_trim_extract_clams.nc'
        rho = []
        icnc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = traj['RHO'].values.flatten()
            rho.extend( vals )
            vals = ( traj['ICN_hom'] + traj['ICN_het'] + traj['ICN_pre'] ).values.flatten()
            icnc.extend( vals )
            vals = traj['T'].values.flatten()
            t.extend( vals )
            vals = traj['RHI'].values.flatten()
            rhi.extend( vals )
    else:
        start_str = 'traj_tst'
        end_str = '_trim_extract.nc'
        rho = []
        icnc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = traj['rho'].values.flatten()
            rho.extend( vals )
            vals = traj['qni'].values.flatten()
            icnc.extend( vals )
            vals = traj['t'].values.flatten()
            t.extend( vals )
            p = traj['p'].values.flatten()
            qv = traj['qv'].values.flatten()
            vals = calc_RHi( vals, p, qv )
            rhi.extend( vals )
    rho = np.array( rho )
    icnc = np.array( icnc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the ICNC value is not 0 or nan
    icnc, t, rhi = filter_icnc_temp( icnc, rho, t, rhi )
    print( icnc.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(icnc)), np.nanmean(np.log10(icnc)), np.nanmax(np.log10(icnc)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, icnc


# Read the ice crystal number concentrations along the CLaMS or ICON trajectories
def read_icnctraj_sub( set_name, time_pt, clams, time0, timef ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    if clams == True:
        start_str = 'cirrus_tst'
        end_str = '_trim_extract_clams.nc'
        rho = []
        icnc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = traj['RHO'].sel( time=slice(time0,timef) ).values.flatten()
            rho.extend( vals )
            vals = ( traj['ICN_hom'] + traj['ICN_het'] + traj['ICN_pre'] ).sel( time=slice(time0,timef) ).values.flatten()
            icnc.extend( vals )
            vals = traj['T'].sel( time=slice(time0,timef) ).values.flatten()
            t.extend( vals )
            vals = traj['RHI'].sel( time=slice(time0,timef) ).values.flatten()
            rhi.extend( vals )
    else:
        start_str = 'traj_tst'
        end_str = '_trim_extract_dt.nc'
        rho = []
        icnc = []
        t = []
        rhi = []
        for file_id in np.arange(1,27):
            filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
            traj = xr.open_dataset( basedir + filename )
            vals = traj['rho'].sel( time=slice(time0,timef) ).values.flatten()
            rho.extend( vals )
            vals = traj['qni'].sel( time=slice(time0,timef) ).values.flatten()
            icnc.extend( vals )
            vals = traj['t'].sel( time=slice(time0,timef) ).values.flatten()
            t.extend( vals )
            p = traj['p'].sel( time=slice(time0,timef) ).values.flatten()
            qv = traj['qv'].sel( time=slice(time0,timef) ).values.flatten()
            vals = calc_RHi( vals, p, qv )
            rhi.extend( vals )
    rho = np.array( rho )
    icnc = np.array( icnc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the ICNC value is not 0 or nan
    icnc, t, rhi = filter_icnc_temp( icnc, rho, t, rhi )
    print( icnc.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(icnc)), np.nanmean(np.log10(icnc)), np.nanmax(np.log10(icnc)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, icnc


# Read the ice water content along the full set of ICON trajectories.
# Below it is assumed that the ice water content has units of kg kg-1. Factor of 10**6 to convert to ppmv
def read_icnctraj_all( set_name, time_pt ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    start_str = 'traj_tst'
    end_str = '.nc'
    rho = []
    icnc = []
    t = []
    rhi = []
    for file_id in np.arange(1,6):
        filename = start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str
        traj = xr.open_dataset( basedir + filename )
        vals = traj['rho'].values.flatten()
        rho.extend( vals )
        vals = traj['qni'].values.flatten()
        icnc.extend( vals )
        vals = traj['t'].values.flatten()
        t.extend( vals )
        p = traj['p'].values.flatten()
        qv = traj['qv'].values.flatten()
        vals = calc_RHi( vals, p, qv )
        rhi.extend( vals )
    rho = np.array( rho )
    icnc = np.array( icnc )
    t = np.array( t )
    rhi = np.array( rhi )

    # Filter out instances where the ICNC value is not 0 or nan
    icnc, t, rhi = filter_icnc_temp( icnc, rho, t, rhi )
    print( icnc.shape, t.shape, rhi.shape )
    print( np.nanmin(np.log10(icnc)), np.nanmean(np.log10(icnc)), np.nanmax(np.log10(icnc)) )
    print( np.nanmin(rhi), np.nanmean(rhi), np.nanmax(rhi) )
    return t, rhi, icnc


# Convert the times in the ICON trajectory files to datetimes
def time_to_datetime( set_name, time_pt, clams ):
    basedir = '/work/bb1018/b380873/traj_output/' + set_name + '/'
    start_str = 'traj_tst'
    end_str = '_trim_extract.nc'
    base_time = datetime( 2017, 8, 6, 0, 0 )
    traj = xr.open_dataset( basedir + start_str + general_prefix(time_pt,8) + '_p001' + end_str )
    traj_times = traj['rtime']
    real_times = np.array( [ base_time + timedelta( seconds=int(t.values) ) for t in traj_times ] )
    for file_id in np.arange(1,31):
        print( file_id )
        traj = xr.open_dataset( basedir + start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + end_str )
        traj['time'] = real_times
        traj.to_netcdf( basedir + start_str + general_prefix(time_pt,8) + '_p' + traj_prefix(file_id) + str(file_id) + '_trim_extract_dt.nc' )

# Temperature dependence of ice water content in both ppmv and g m-3 from
# Martina Kraemer in-situ climatology (see ACP 2020 Cirrus Guide)
def martina_T_IWC_line():
    mh2o = 18 # Molar mass of water [g mol-1]
    R = 8.31 # Gas constant [J K-1 mol-1]
    rho_ice = 0.92 # Density of ice [g cm-3]

    T_min = 180
    T_max = 240
    Tgrid = np.arange(T_min, T_max+1)

    p0 = 0.000289319
    p1 = 2.66906
    p2 = -247.724
    pmean = p0*Tgrid**p1 + p2

    a0 = -54.3215
    a1 = 0.986506
    a2 = 2.87180 - 0.3
    IWC_ppm_min = 10**(a0*a1**Tgrid + a2)
    factor = (pmean*10**2 / Tgrid) * mh2o / R
    IWC_gm3_min = IWC_ppm_min*10**(-6)*factor

    a0 = -396248
    a1 = 0.937089
    a2 = 2.04491 + 0.4
    IWC_ppm_max = 10**(a0*a1**Tgrid + a2)
    IWC_gm3_max = IWC_ppm_max*10**(-6)*factor

    a0 = -68.4202
    a1 = 0.983917
    a2 = 2.81795
    IWC_ppm_med = 10**(a0*a1**Tgrid + a2)
    IWC_gm3_med = IWC_ppm_med*10**(-6)*factor

    return Tgrid, IWC_ppm_min, IWC_ppm_max, IWC_ppm_med

# Percentiles of temperature and ice crystal number from Martina Kraemer
# in-situ climatology (see ACP 2020 Cirrus Guide)
def martina_T_Ni_perc():
    T_perc = [180.0, 181.0, 182.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 189.0,
              190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 199.0,
              200.0, 201.0, 202.0, 203.0, 204.0, 205.0, 206.0, 207.0, 208.0, 209.0,
              210.0, 211.0, 212.0, 213.0, 214.0, 215.0, 216.0, 217.0, 218.0, 219.0,
              220.0, 221.0, 222.0, 223.0, 224.0, 225.0, 226.0, 227.0, 228.0, 229.0,
              230.0, 231.0, 232.0, 233.0, 234.0, 235.0, 236.0, 237.0, 238.0, 239.0,
              240.0, 241.0, 242.0, 243.0, 244.0]

    Ni_perc_10 = [0.00000,  0.00000,  0.00000,  0.00000,  0.00000, 0.000688500, 0.000953300,
     0.000936300, 0.000922000, 0.000847600, 0.000953500, 0.000730000, 0.000762400,  0.001146,
     0.00108980,  0.00129870,  0.00126470,  0.00121120,  0.00152730,  0.00146800,  0.00140500,
     0.00113200,  0.00127900,  0.00110685,  0.00115230,  0.00119000,  0.00225140,  0.00325240,
     0.00167000,  0.00217000,  0.00220000,  0.00542000,  0.00479000,  0.00211000,  0.00322000,
     0.00389000,  0.00361000,  0.00384000,  0.00566685,  0.00286000,  0.00349000,  0.00212786,
     0.00230000,  0.00171000,  0.00131000,  0.00278000,  0.00289000,  0.00269000,  0.00438000,
     0.000440000,  0.00125000,  0.00542000, 0.000690000,  0.00166000, 0.000290000,  0.00202000,
     0.000730000, 0.000680000, 0.000290000,  0.00122000,  0.00214000, 0.000830000, 0.000350000,
     0.000390000,  0.00167000]

    Ni_perc_25 = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000, 0.00245500, 0.0172171, 0.0122320,\
     0.00413700,  0.00280700,  0.00324000,  0.00236080,  0.00282780,  0.00401558,  0.00333000, 0.00492760,\
     0.00438970,  0.00467810,  0.005057,  0.00532560,  0.00451000,  0.00348760,  0.00470770,  0.00473030,\
     0.00383840,  0.00356480, 0.0134855, 0.0136300, 0.0102500, 0.0140300, 0.0100200, 0.0227800, 0.0169100,\
     0.00751000, 0.0145300, 0.0124700, 0.0108000, 0.0110100, 0.0138500,  0.00731000,  0.00882000, 0.00817000,\
     0.00769000, 0.00706000, 0.00582000, 0.0110200, 0.00971000, 0.00639000, 0.0121600, 0.00475000, 0.00593000,\
     0.0144500,  0.00294000,  0.00698000,  0.00214000,  0.00780000,  0.00465000,  0.00447000,  0.00192000,\
     0.00471000,  0.00940000,  0.00966000,  0.00203000,  0.00376000, 0.0116700]

    Ni_perc_50 = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  0.0342310,  0.0554229,  0.0555400,  0.0448850,
     0.03239, 0.0300830, 0.0218721, 0.0173452, 0.0181934, 0.0181052, 0.0236230, 0.0183167, 0.0189614, 0.0178192,
     0.0209973,  0.0140938,  0.0165819,  0.0310280,  0.0356804,  0.0156300,  0.0155887,  0.0464656,  0.0400700,
     0.0410100,  0.0576500,  0.0407700,  0.0583300,  0.0376800,  0.0296700,  0.0480400,  0.0371500,  0.0483200,
     0.0353600,  0.0471700,  0.0304100,  0.0300200,  0.0358100,  0.0231900,  0.0351600,  0.0496600,  0.0585400,
     0.0397084,  0.0287600,  0.0467500,  0.0226000,  0.0297900,  0.0510800,  0.0115800,  0.0355606, 0.00916000,
     0.0224800,  0.0123600,  0.0142200,  0.0187400,  0.0218200,  0.0309900,  0.0685400,  0.0171900,  0.0219906,
     0.0514000 ]

    Ni_perc_75 = [0.00000, 0.00000, 0.00000, 0.00000, 0.00000,  0.0815330,  0.0968640, 0.131302, 0.126379, 0.110781,
     0.0955200,  0.0746895,  0.0560400,  0.0639110,  0.0740875,  0.0638451,  0.0536802,  0.0524256,  0.0497950,
     0.0551680,  0.0359832,  0.0531353,  0.0935090,0.108604,  0.0564388,  0.0694435, 0.128646, 0.129450, 0.134870,
     0.185190, 0.142080, 0.188950, 0.125880, 0.0812500, 0.122560, 0.0762200, 0.135370, 0.0998100, 0.192460, 0.147160,
     0.134050, 0.143540, 0.0833600, 0.102980, 0.220880, 0.244164, 0.122840, 0.123080, 0.112860, 0.115380, 0.131210,
     0.171730, 0.0499000,  0.109830, 0.0485000,  0.100990, 0.0528800, 0.0737600,  0.122814, 0.0806000,  0.148105,
     0.216730,  0.184730,  0.133900,  0.151260 ]

    Ni_perc_90 = [0.00000,   0.00000,   0.00000,   0.00000,   0.00000,  0.188489,  0.149230,  0.399597,  0.405776,
     0.370650,  0.207402,  0.169595,  0.146242,  0.148380,  0.222482,  0.181361,  0.132376,  0.132272,  0.150940,
     0.143223, 0.0920229,  0.130170,  0.240383,  0.236314,  0.145874,  0.236818,  0.258815,  0.298060,  0.326656,
     0.467660,  0.384610,  0.855060,  0.418030,  0.245690,  0.278820,  0.204190,  0.273950,  0.279050,  0.704930,
     0.503160,  0.473830,  0.336880,  0.200810,  0.239930,  0.489860,  0.580320,  0.332100,  0.317850,  0.226790,
     0.365330,  0.623790,  0.446540,  0.133150,  0.209736,  0.131680,  0.408590,  0.133560,  0.301570,  0.920610,
     0.283135,  0.421600,  0.492530,   1.46023,  0.607900,  0.308620]

    return T_perc, Ni_perc_10, Ni_perc_25, Ni_perc_50, Ni_perc_75, Ni_perc_90


# Percentiles of temperature and ice crystal effective radius [micron] from Martina Kraemer
# in-situ climatology (see ACP 2020 Cirrus Guide)
def martina_T_ri_perc():
    T_perc = [180.0, 181.0, 182.0, 183.0, 184.0, 185.0, 186.0, 187.0, 188.0, 189.0,
              190.0, 191.0, 192.0, 193.0, 194.0, 195.0, 196.0, 197.0, 198.0, 199.0,
              200.0, 201.0, 202.0, 203.0, 204.0, 205.0, 206.0, 207.0, 208.0, 209.0,
              210.0, 211.0, 212.0, 213.0, 214.0, 215.0, 216.0, 217.0, 218.0, 219.0,
              220.0, 221.0, 222.0, 223.0, 224.0, 225.0, 226.0, 227.0, 228.0, 229.0,
              230.0, 231.0, 232.0, 233.0, 234.0, 235.0, 236.0, 237.0, 238.0, 239.0,
              240.0, 241.0, 242.0, 243.0, 244.0]

    Ri_perc_10 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 4.07087, 4.39798, 3.55754, 3.23484, 3.94133, 4.39245,
                  5.30966, 6.07553, 6.63280, 7.02325, 7.62439, 8.07797, 8.64927, 8.89180, 9.24180, 11.9808, 12.3634,
                  10.2635, 12.6437, 12.79, 12.1359, 11.2031, 11.1814, 9.20913, 8.75288, 8.97744, 7.70370, 13.0291,
                  14.3112, 14.1418, 15.3197, 13.3327, 16.8432, 16.2016, 16.9101, 19.6084, 19.7614, 21.7810, 21.6510,
                  20.4422, 21.4036, 25.5874, 24.6495, 24.5772, 22.0104, 22.2741, 28.5583, 24.8485, 29.5459, 24.0526,
                  25.5810, 22.5049, 25.7283, 27.5581, 24.0726, 26.9062, 23.5407, 26.2579, 25.9905, 22.8270]

    Ri_perc_25 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 5.00621, 5.83743, 4.72092, 4.92424, 5.90791, 6.71322, 8.25986,
                  10.2003, 9.78089, 10.1030, 12.16, 12.15, 13.2381, 13.4457, 13.9009, 15.4059, 16.3154, 15.5254, 16.2726,
                  16.6956, 16.8693, 14.7705, 17.7920, 14.2711, 13.0070, 13.9881, 13.7110, 15.9943, 16.9833, 17.7571, 18.4625,
                  17.3808, 22.3364, 20.9339, 21.3660, 24.45, 25.6418, 28.1055, 27.9766, 26.0410, 26.0719, 32.4834, 29.5486,
                  31.7273, 28.3131, 28.5239, 35.1795, 31.6518, 37.6468, 30.8871, 32.1236, 33.1780, 35.5972, 32.6109, 33.2352,
                  33.9963, 29.7821, 30.1006, 33.1120, 34.3615]

    Ri_perc_50 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 7.18534, 7.34107, 6.80043, 7.60545, 9.78124, 11.9980, 15.3632,
                  16.6502, 16.1359, 16.5496, 17.7878, 16.9082, 18.1071, 18.9765, 19.0335, 18.4360, 19.9817, 19.7785, 20.4125,
                  21.7341, 22.4894, 19.5016, 22.2969, 18.8649, 18.2578, 17.8652, 18.3978, 19.4475, 20.5430, 22.0177, 21.8979,
                  24.2725, 29.4083, 26.8723, 28.7156, 30.7340, 33.0439, 37.1119, 38.6392, 32.8211, 33.2526, 42.2811, 41.4803,
                  42.3508, 37.4024, 35.9896, 41.4169, 42.7797, 47.6127, 42.6400, 41.4618, 47.8127, 49.4781, 43.6934, 45.8360,
                  43.6462, 44.1776, 42.5178, 45.0298, 44.6836]

    Ri_perc_75 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 17.5850, 11.2590, 15.2311, 17.0926, 18.2015, 18.6361, 20.1588,
                  21.3517, 20.7865, 21.4096, 22.0992, 21.1252, 23.6307, 24.2060, 23.7612, 22.5575, 24.6358, 24.8467, 25.4546,
                  27.1219, 28.5598, 25.5906, 28.9455, 25.9011, 24.7818, 23.7859, 22.7497, 23.4512, 24.5414, 28.9339, 27.3822,
                  31.9649, 39.6358, 35.1041, 40.2067, 40.3732, 43.2939, 49.6404, 54.7999, 43.1115, 43.0584, 53.5447, 59.1823,
                  56.5753, 47.5540, 47.8308, 54.1276, 56.5386, 65.2121, 56.7153, 57.1806, 67.0175, 63.2496, 57.1532, 61.6903,
                  56.1368, 58.4143, 61.6989, 61.9442, 59.3760]

    Ri_perc_90 = [0.0000, 0.0000, 0.0000, 0.0000, 0.0000, 23.5772, 22.7807, 20.8354, 21.3845, 22.5355, 22.8195, 25.0678,
                  26.8571, 27.2934, 29.2247, 25.9860, 26.3874, 29.3436, 29.3004, 28.6535, 27.0905, 29.0920, 29.8917, 31.2387,
                  32.6968, 34.4625, 32.4538, 35.2411, 35.5832, 31.6078, 30.2126, 27.5559, 27.8971, 30.1043, 46.0018, 46.0378,
                  40.6349, 50.1515, 47.9232, 52.0012, 50.7079, 59.7474, 65.0421, 68.0434, 60.1175, 60.2523, 64.0725, 69.2324,
                  69.9230, 61.6686, 61.8274, 65.6581, 67.6866, 74.6206, 70.0884, 71.7930, 77.6144, 73.6540, 71.3246, 77.6974,
                  69.7655, 76.7911, 78.3845, 76.3504, 73.3093]

    return T_perc, Ri_perc_10, Ri_perc_25, Ri_perc_50, Ri_perc_75, Ri_perc_90
