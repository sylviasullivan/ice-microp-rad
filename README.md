# ICON_vis
This repository contains a variety of python scripts to visualize, regrid, and analyze ICON model output from a limited-area domain simulation. Scripts are grouped into those that visualize the domain, dynamical fields, ice / cloud hydrometeor fields, cloud optics values, radiative values, and trajectories. There are also remapping and utility directories.

*fig2_reproduced.py* - reproduce Figure 2 of Lee et al. ACP doi: 10.5194/acp-19-11803-2019, i.e. qv, qi, temperature, potential temperature, and zonal wind profiles from StratoClim Flight 7 versus the model

CloudSat/ - scripts to visualize the 2CICE ice water content (IWC) product, derived from CloudSat overpasses, with the ICON output IWC
- *2CICE_ZL.py* - generate altitude-lat/lon swathes of 2CICE and ICON IWC 
- *CloudSat_read.py* - utility to read in the CloudSat data format, modified from the default from Geroge (CloudSat_read_George.py)

domain/ - scripts to visualize the simulation domain and vertical grids
- *flight_track.py* - latitude over time of Flight 7
- *topography.py* - map of simulation domain, Flight 7 track, two CloudSat overpasses
- *vgrid_spacings.py* - layer thickness versus altitude for the different vertical grids
- *vgrid_vis.py* - level count versus altitude for the increased and default vertical grids

dynamics/ - script to visualize the flow fields
- *MMCR_ICON.py* - comparison of the millimeter cloud radar and ICON vertical velocities
- *wComp.py* - spatial comparison of ERA5 and ICON-simulated pressure velocities at 250 and 500 hPa
- *wDist_subdomains.py* - visualization of data from Silvia Bucci
- *winds.py, wind_lowlevel.py, winds_hilevel.py* - visualization of wind fields at three levels (850, 500, 250 hPa) or just one level
- *wT-profiles.py* - deviations in daily-mean, domain-mean pressure velocity and temperature profiles from the 4 simulations

ice/ - scripts to visualize the ice fields
- *DARDAR_Nice.py* - global climatology of DARDAR ice crystal number values
- *diurnal-IWPdist.py, diurnal-IWPdist2.py* - diurnally separated probability distributions of IWP (absolute or difference)
- *IWPdist.py* - log-linear probability distribution of IWP values for 1-mom and 2-mom simulations
- *IWPdist_land-ocean.py* - probability distribution of IWP over land versus ocean
- *IWPdist_olr.py* - probability distributions of IWP filtered for different OLR bins
- *IWPdist_precip.py* - probability distributions of IWP filtered for different precipitation accumulation bins
- *IWP-wDist.py* - probability distributions of IWP, w_250hPa, w_500hPa values for 1-mom and 2-mom simulations
- *qi-field.py* - spatial map of ice mass mixing ratio
- *q\*-profiles.py* - daily-mean domain-mean profiles of ice, snow, graupel, or cloud droplet mass mixing ratio for the 1-mom and 2-mom simulations
- *qi_video.py* - make a video of the ice mass mixing ratio evolution over time
- *tqi_video.py* - make a video of the column-integraed ice content over time

optics/ - scripts to visualize the cloud optics fields
- *cldopt_vis.py* - visualize the 2D surface plot of mass absorption coefficient (or single-scattering albedo, asymmetry parameter with modifications) as a function of ice crystal effective radius and wavelength
- *cldopt_vis2.py* - visualize the mass absorption coefficient, single-scattering albedo, and asymmetry parameters as a function of wavelength, colored for different effective radii
- *cldopt_vis3.py* - Baum-Yang optical properties visualization

rad/ - scripts to visualize the radiative fields
- *OLR_v_IWP.py* - visualizes the TOA (outgoing) longwave flux as a function of ice water path (IWP), needs modification re: files read
- *flux-profiles.py* - domain-mean, daily-mean longwave and shortwave profiles in different simulations
- *heating-profiles.py* - domain mean, daily-mean longwave and shortwave heating profiles in different simulations 
- *clrsky-profiles.py* - domain mean, daily mean clear-sky longwave and shortwave heating profiles in different simulations
- *diurnal-heating.py* - heating profiles for each hour throughout the day
- *olrComp.py, olrComp2.py* - spatial comparisons of outgoing longwave radiation (OLR) (between CERES, ERA5, and ICON sims *or* between the 4 ICON sims)
- *olrDist.py* - probability density comparisons of outgoing longwave radiation (OLR) between CERES, ERA5, and ICON, shows these individually with their KL div
- *olrDist2.py, olrDist3.py* - probability density comparisons of outgoing longwave radiation (OLR) (between CERES, ERA5, and ICON sims *or* between the 4 ICON sims)
- *olrDist_IWP-w.py* - OLR probability density filtered by bins of IWP, w_250hPa, or w_500hPa
- *olrDist_w.py* - 
- *olr_video3.py, olr_video_CERES.py, olr_video_ERA5.py* - make a video of the OLR field from the ICON simulation, CERES measurements, or ERA5 values over time

remapping/ - scripts to perform vertical or horizontal remapping
- *AP_to_HL.sh* - uses the ap2hl command to vertically regrid files to altitudes from model levels, given file prefixes and first and last timesteps (part2 variable can be used to extract certain variables as well)
- *AP_to_PL.sh, AP_to_PL2.sh* - uses the ap2pl command to vertically regrid files to pressures from model levels, given file prefixes and first and last timesteps (part2 variable can be used to extract certain variables as well)
- *meanPL.py* - generate a list of the mean pressures per model level and save this to a text file (that can later be used by *AP_to_PL\*.sh*)
- *submit_ap2pl.sh* - slurm version of *AP_to_PL\*.sh*
- *submit_gencon.sh, submit_gendis.sh* - slurm script to generate variable-conserving or distance-weighted remapping weights
- *submit_griddes.sh* - slurm script to generate a grid description text file
- *submit_remapcon.sh, submit_remapdis.sh* - slurm script to perform variable-conserving or distance-weighted remapping 

traj/ - scripts to visualize the trajectory outputs
- *lagranto_vis.py* - attempt to use the lagranto.plotting package
- *qv_time.py* - compare the mean specific humidity along CLaMS-ice and ICON trajectories
- *traj_nc.py* - trajectory file postprocessing script
- *traj_time.py, traj_time2.py* - visualize the ice mass mixing ratio in the CLaMS-ice versus ICON trajectory over time (the first also shows filtering for same trajectory mean temperature)
- *traj_psd.py* - some utilities from Christian Rolf to smooth and calculate power spectral densities from trajectory output
- *traj_vis.py* - visualize the trajectories colored by altitude over a land type-colored map
- *psd_driver.py* - visualize the temperature and vertical velocity power spectral densities from trajectory output

utilities 
- *convertTXTNC.py* - generate nc files from the grid description txt file
- *KL_div.py* - calculate the Kullback-Leibler divergence between two distributions
- *z_from_ml.py* - calculate the approximate altitude for a given model level at lat/lon using the external parameter topography
