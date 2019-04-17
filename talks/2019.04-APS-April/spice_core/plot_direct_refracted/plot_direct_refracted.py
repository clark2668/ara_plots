import matplotlib.pyplot as plt
import numpy as np
from numpy import genfromtxt
import time
from NuRadioMC.SignalProp import analyticraytraycing as ray
from NuRadioMC.utilities import units, medium
import logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('raytracing')
from pylab import setp
import matplotlib.patches as patches

fig = plt.figure(figsize=(2*11.5,8.5))
sizer=32

####
## Ray Tracing
###
fig.subplots_adjust(wspace=2.5)
ax = fig.add_subplot(1,2,1)

x1 = [-2382 * units.m, -1130. * units.m]  # pulser position
x2 = [0., -200. * units.m]  # ARA antanna
r = ray.ray_tracing_2D(medium.southpole_simple())
solution = r.find_solutions(x1, x2)
dirY, dirZ = r.get_path(x1, x2, solution[0]['C0'])
refY, refZ = r.get_path(x1, x2, solution[1]['C0'])

ax.plot(-dirY/ units.km, dirZ/ units.km,'-',linewidth=5.0,color='royalblue',label=r'Direct Solution')
ax.plot(-refY/ units.km, refZ/ units.km,'-',linewidth=5.0,color='firebrick',label=r'Refracted Solution')

ax.plot(0,-0.2,"o",color='black',markersize=20)
ax.plot(2.382,-1.130,"D",color='black',markersize=20)

vert_line_x=np.array([2.382,2.382])
vert_line_y=np.array([-1.130,0])
ax.plot(vert_line_x,vert_line_y,'--',linewidth=4.0,color='black',label=r'hole')


# ax.arrow( 2.598, -1.110, -0.167, 0.06,  head_width=0.07, head_length=0.1, fc='royalblue', ec='royalblue')
# ax.arrow( 2.601, -0.9957, -0.153, 0.0747,  head_width=0.07, head_length=0.1, fc='firebrick', ec='firebrick')


ax.set_xlim(-0.1, 2.5)
ax.set_ylim(-1.3, 0.0)
ax.set_xticks([0,0.5,1,1.5,2,2.5])
# ax.set_yticks([-1.4,-1.,-.6,-.2])

ax.tick_params(labelsize=sizer-3)
ax.set_xlabel('Horizontal Distance (km)', size=sizer)
ax.set_ylabel('Depth (km)', size=sizer)

ax.text(0.42, 0.33, 'Direct\nSolution',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='royalblue', fontsize=29)

ax.text(0.62, 0.6, 'Refracted\nSolution',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='firebrick', fontsize=29)

ax.text(0.01, 0.70, 'ARA2',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=29)

ax.text(0.72, 0.04, 'SPIce Pulser',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=29)

###
# Waveform
##

ax2 = fig.add_subplot(1,2,2)

data = np.genfromtxt("waveform.txt",delimiter=',',skip_header=1,names=['time','volts'])
split=1200
ax2.plot(data['time'][:split],data['volts'][:split]/1000.,'-',linewidth=2.0,color='royalblue',label=r'Direct Solution')
ax2.plot(data['time'][split:],data['volts'][split:]/1000.,'-',linewidth=2.0,color='firebrick',label=r'Refracted Solution')

ax2.set_xlim(-50, 700)
ax2.set_ylim(-1, 1)
# ax2.set_xticks([100,200,300,400,500])

ax2.tick_params(labelsize=sizer-3)
ax2.set_xlabel('Time (ns)', size=sizer)
ax2.set_ylabel('Voltage (V)', size=sizer)


ax2.text(0.15, 0.8, 'Direct\nPulse',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax2.transAxes,
        color='royalblue', fontsize=29)

ax2.text(0.75, 0.8, 'Refracted\nPulse',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax2.transAxes,
        color='firebrick', fontsize=29)

fig.tight_layout()
fig.savefig("spice_to_a2.png",dpi=300)

