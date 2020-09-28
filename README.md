# ICON_vis

CloudSat/ - scripts to visualize the 2CICE ice water content (IWC) product, derived from CloudSat overpasses, with the ICON output IWC
- 2CICE_ZL.py - generate altitude-lat/lon swathes of 2CICE and ICON IWC 
- CloudSat_read.py - utility to read in the CloudSat data format, modified from the default from Geroge (CloudSat_read_George.py)

domain/ - scripts to visualize the simulation domain and vertical grids
- flight_track.py - latitude over time of Flight 7
- topography.py - map of simulation domain, Flight 7 track, two CloudSat overpasses
- vgrid_spacings.py - layer thickness versus altitude for the different vertical grids
- vgrid_vis.py - level count versus altitude for the increased and default vertical grids

