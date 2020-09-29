import cartopy.crs as ccrs
import cartopy.feature as cfeature
#from lagranto.plotting import plot_trajs
import matplotlib.pyplot as plt

crs = ccrs.Stereographic(central_longitude=65,central_latitude=45)

fig = plt.figure()
ax = plt.axes(projection=crs)
land_50m = cfeature.NaturalEarthFeature('cultural', 'admin_0_countries',
                    '50m', edgecolor='gray', facecolor='none', linewidth=0)
ax.add_feature(land_50m)
ax.set_extent([50,115,-5,40])
plt.show()
