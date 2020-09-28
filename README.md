# ICON_vis

CloudSat/ - scripts to visualize the 2CICE ice water content (IWC) product, derived from CloudSat overpasses, with the ICON output IWC
- *2CICE_ZL.py* - generate altitude-lat/lon swathes of 2CICE and ICON IWC 
- *CloudSat_read.py* - utility to read in the CloudSat data format, modified from the default from Geroge (CloudSat_read_George.py)

domain/ - scripts to visualize the simulation domain and vertical grids
- *flight_track.py* - latitude over time of Flight 7
- *topography.py* - map of simulation domain, Flight 7 track, two CloudSat overpasses
- *vgrid_spacings.py* - layer thickness versus altitude for the different vertical grids
- *vgrid_vis.py* - level count versus altitude for the increased and default vertical grids

dynamics/ - script to visualize the flow fields
- **

ice/ - scripts to visualize the ice fields

rad/ - scripts to visualize the radiative fields
- *OLR_v_IWP.py* - visualizes the TOA (outgoing) longwave flux as a function of ice water path (IWP), needs modification re: files read
- *flux-profiles.py* - domain-mean, daily-mean longwave and shortwave profiles in different simulations
- *heating-profiles.py* - domain mean, daily-mean longwave and shortwave heating profiles in different simulations 
- *diurnal-heating.py* - heating profiles for each hour throughout the day
- *olrComp.py, olrComp2.py* - spatial comparisons of outgoing longwave radiation (OLR) (between CERES, ERA5, and ICON sims *or* between the 4 ICON sims)
- *olrDist.py* - probability density comparisons of outgoing longwave radiation (OLR) between CERES, ERA5, and ICON, shows these individually with their KL div
- *olrDist2.py, olrDist3.py* - probability density comparisons of outgoing longwave radiation (OLR) (between CERES, ERA5, and ICON sims *or* between the 4 ICON sims)
- *olrDist_IWP-w.py* - OLR probability density filtered by bins of IWP, w_250hPa, or w_500hPa
- *olrDist_w.py* - 
- *olr_video3.py, olr_video_CERES.py, olr_video_ERA5.py* - make a video of the OLR field from the ICON simulation, CERES measurements, or ERA5 values over time

utilities 
- *KL_div.py* - calculate the Kullback-Leibler divergence between two distributions
