# -*- coding: utf-8 -*-
import numpy as np #import numpy
import matplotlib.pyplot as plt
import math
import util

mode = 'km'

icecube_center = util.get_thing('I3Center', mode)
print("IceCube Center {}".format(icecube_center))

xs, ys = util.get_icecube_string_locations(mode)
xs+=icecube_center[0]
ys+=icecube_center[1]

SPIce = util.get_thing('SPIce', mode)
print("SPIce {}".format(SPIce))

IC22 = util.get_thing('IC22', mode)
print("IC22 {}".format(IC22))

IC1 = util.get_thing('IC1', mode)
print("IC1 {}".format(IC1))

RTP = util.get_thing('RTP', mode)
print("RTP {}".format(RTP))

A2 = util.get_thing('A2', mode)
print("A2 {}".format(A2))

fig = plt.figure(figsize=(5*1.5,5*1.5))
ax=fig.add_subplot(111)

ax.plot(xs,ys,'o',color='gray',linewidth=0, markersize=5, alpha=0.2)
ax.plot(A2[0],A2[1], "o", label="A2")
ax.plot(SPIce[0],SPIce[1],"kx", label="SPIce")
ax.plot(IC22[0],IC22[1],"k+", label="IC22")
ax.plot(IC1[0],IC1[1],"kv", label="IC1")
ax.plot(RTP[0], RTP[1], 'k^', label='RTP')
	
ax.set_xlabel('Easting (km)')
ax.set_ylabel('Northing (km)')

ax.set_xlim(10,15)
ax.set_ylim(12,17)
ax.grid()
ax.legend(loc='lower right')

fig.savefig("a2_spice_icecube.png",edgecolor='none',bbox_inches="tight")