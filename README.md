# ice-microp-rad
This repository contains a variety of Python scripts to visualize, regrid, and analyze ICON v 2.3.0 model output from limited-area domain and trajectory simulations. The project focuses on the interaction of ice microphysics with radiation. Scripts are grouped into those that visualize the domain, dynamical fields, ice / cloud hydrometeor fields, cloud optics values, radiative values, and trajectories. There are also remapping and utility directories.

*Figures_article1.ipynb*, *SI_article1.ipynb* - iPython notebooks that pull from the functions below to generate the 4 main and 5 supplemental figures of the [Communications Earth & Environment publication]{https://www.nature.com/articles/s43247-021-00206-7}

*Figures_article2.ipynb*, *SI_article2.ipynb* - iPython notebooks that pull from the functions below to generate figures for a publication on the trajectories in ICON and CLaMS

*Postprocessing1_article2.ipynb*, *Postprocessing2_article2.ipynb* - iPython notebooks that calculate probability distributions from the CLaMS and ICON trajectory output. Postprocessing1 combines all trajectory files and calculates PDFs of cloud ice mass and number, also versus temperature and RHi. Postprocessing2 digitizes and saves ICON-CLaMS differences in cloud ice mass and number along either an absolute or normalized time axis.

*Analysis1_article2.ipynb*, *Analysis2_article2.ipynb* - iPython notebooks that act as workspaces to refine figures prior to Figures_article2. Analysis1_article2 works with StratoClim measurements and ICON synthetic flight tracks. Analysis2_article2 works with differences between ICON and CLaMS trajectory output.

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
- *reff_fix.ipynb* - visualize mean profiles of droplet and crystal effective radius or number concentrations

remapping/ - scripts to perform vertical or horizontal remapping
- *AP_to_HL.sh* - uses the ap2hl command to vertically regrid files to altitudes from model levels, given file prefixes and first and last timesteps (part2 variable can be used to extract certain variables as well)
- *AP_to_PL.sh, AP_to_PL2.sh* - uses the ap2pl command to vertically regrid files to pressures from model levels, given file prefixes and first and last timesteps (part2 variable can be used to extract certain variables as well)
- *meanPL.py* - generate a list of the mean pressures per model level and save this to a text file (that can later be used by *AP_to_PL\*.sh*)
- *submit_ap2pl.sh* - slurm version of *AP_to_PL\*.sh*
- *submit_gencon.sh, submit_gendis.sh* - slurm script to generate variable-conserving or distance-weighted remapping weights
- *submit_griddes.sh* - slurm script to generate a grid description text file
- *submit_remapcon.sh, submit_remapdis.sh* - slurm script to perform variable-conserving or distance-weighted remapping 

traj/ - scripts to visualize the trajectory outputs
- *CLaMS-negRHi-exploration.ipynb* - visualizes time series of CLaMS RHi versus various cloud variables along the ICON or CLaMS trajectories to understand the negative RHi values produced when using the SH_mix tendency
- *diff_abs_lifetime.py*, *diff_norm_lifetime.py* - calculate and save the ICON-CLaMS difference in qi or Ni versus a digitized absolute or normalized time
- *flighttrack.py* - three functions related to StratoClim data: (1) read in values, (2) bin them altitudinally, (3) retain only whole-second data
- *icetraj.py* - functions in a preliminary attempt to extract and filter qi, Ni, T, and RHi along trajectories (replaced by Postprocessing1_article2.ipynb); includes visualization of the T-IWC climatology and T-Ni-ri percentiles from Martina Krämer
- *psd_driver.py* - visualize the temperature and vertical velocity power spectral densities from trajectory output
- *lagranto_vis.py* - attempt to use the lagranto.plotting package
- *t_pdf.py* - PDF of temperature and updraft fluctuations along the ICON trajectories
- *total_water.ipynb* - one-to-one or relative difference plots of total water (qv+qi) in ICON versus CLaMS-ice simulations
- *traj_nc.py* - trajectory file postprocessing script
- *traj_psd.py* - some utilities from Christian Rolf to smooth and calculate power spectral densities from trajectory output
- *traj_vis.py* - visualize the trajectories colored by altitude over a land type-colored map

traj/syntraj/ - scripts to work with synthetic trajectories from ICON high-time-resolution output
- *collocateme.py* - driver for collocateSim to produce the trajectory with minimized Euclidean distance
- *collocateSim.py* - find the ICON values with minimized Euclidean distance from the in-situ measurement
- *extractme.py* - driver for extractSim to produce a series of synthetic trajectories in a spatiotemporal "bubble" around the in-situ measurement
- *extractSim.py* - find a random subset of ICON values within a spatiotemporal "bubble" around the in-situ measurement to constitute synthetic trajectories
- *pinpointme.py* - driver for pinpointSim to produce the trajectory with values closest to those measured along the flight track
- *pinpointSim.py* - find the ICON values closest to those measured along the flight track
- *statme_\*.py* - driver for the corresponding *syn_traj_stats_\*.npy* to calculate statistics in altitude bins along synthetic trajectories
- *syntraj_stats_\*.py* - calculate statistics in altitude bins along the synthetic trajectories (*fixed* with the same element number per bin as the measurements)

utilities/  
- *convertTXTNC.py* - generate nc files from the grid description txt file
- *densityPlot.py* - two functions to automatize plotting of 2D density for qi/Ni versus T/RHi and store means and medians of these variables
- *KL_div.py* - calculate the Kullback-Leibler divergence between two distributions
- *plotting_utilities.py* - series of functions to make nice axes, generate a color dictionary, open certain files, etc.
- *thermodynamic_functions.py* - calculate saturation vapor pressure, potential temperature, RH, etc.
- *z_from_ml.py* - calculate the approximate altitude for a given model level at lat/lon using the external parameter topography
