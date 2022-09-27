import os
import numpy as np
import ROOT
from matplotlib import pyplot as plt
import itertools

ROOT.gSystem.Load(os.environ.get('ARA_UTIL_INSTALL_DIR')+"/lib/libAraEvent.so")
ara_geom = ROOT.AraGeomTool.Instance()

markers = itertools.cycle(('o', 's', '^', 'v', '>', '<'))

fig = plt.figure(figsize=(8,8))
ax = plt.subplot(121, projection='polar')
ax2 = plt.subplot(122)

stations = np.linspace(1,5, 5)
for st in stations:
    st = int(st)
    stationVector = ara_geom.getStationVector(st)
    lon = ara_geom.getLongitudeFromArrayCoords(stationVector[1], stationVector[0], 2011) # get the longitude
    lat = ara_geom.getGeometricLatitudeFromArrayCoords(stationVector[1], stationVector[0], 2011) # get the 'Geometric' latitude. We are going to compare with weather balloon 'XYZ' coordinates
    lat = np.radians(lat)
    lon = np.radians(lon) 
    # print(f'A{st} coord. Lat: {np.degrees(lat)} deg, Lon: {np.degrees(lon)} deg') 
    # print(stationVector[0], stationVector[1], stationVector[2])
    
    x = stationVector[0]
    y = stationVector[1]
    the_mark = next(markers)
    the_label = f'A{st}'
    ms = 10
    ax.plot([lon], [lat], marker=the_mark, label=the_label, markersize=ms)
    ax2.plot([x], [y], marker=the_mark, label=the_label, markersize=ms)

ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)
ax.set_rlabel_position(45)

ax2.set_aspect('equal')
ax2.legend()



fig.tight_layout()
fig.savefig('map.png')
