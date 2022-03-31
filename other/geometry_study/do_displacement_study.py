# -*- coding: utf-8 -*-
from re import A
from shutil import which
from sys import displayhook
import numpy as np #import numpy
import matplotlib.pyplot as plt
import math
import util
import copy

mode = 'meters'

SPIce = util.get_thing('SPIce', mode)
RTP = util.get_thing('RTP', mode)
A2 = util.get_thing('A2', mode)

diff_RTP = RTP - A2
true_angle_rtp =  util.get_angle(diff_RTP)
print("Nominal RTP Angle {}".format(true_angle_rtp))

diff_spice = SPIce - A2
true_angle_spice =  util.get_angle(diff_spice)
print("Nominal Spice Angle {}".format(true_angle_spice))

# MYL sees the RTP at 257.1 (0->360), which means -103 in -180->180
MYL_A2_RTP_AZI = -102.9
MYL_A2_RTP_AZI = util.transform_map_to_global(MYL_A2_RTP_AZI)
print("MYL Observation {}".format(MYL_A2_RTP_AZI))

# Justin sees the Spice pulser at 255.1 (0->360), which means -104.9
JF_A2_SPICE_AZI = -104.90
JF_A2_SPICE_AZI = util.transform_map_to_global(JF_A2_SPICE_AZI)
print("JF Observation {}".format(JF_A2_SPICE_AZI))

# simple 1D scan
do_1D_scan  = False
if do_1D_scan:

    which_dir = 'horizontal'
    # which_dir = 'vertical'
    sources = ['spice', 'rtp']
    displacements = np.linspace(-1500, 1500, 101)
    angles = {}
    angles['spice'] = []
    angles['rtp'] = []
    for d in displacements:
        A2_temp = copy.copy(A2)
        if which_dir =='vertical':
            A2_temp[1] += d
        elif which_dir =='horizontal':
            A2_temp[0] += d
        diff_rtp = RTP - A2_temp
        diff_spice = SPIce - A2_temp
        angle_rtp = util.get_angle(diff_rtp)
        angle_spice = util.get_angle(diff_spice)
        # print("D {}, Angle {}".format(d, angle))
        angles['rtp'].append(angle_rtp)
        angles['spice'].append(angle_spice)

    fig = plt.figure(figsize=(5, 5))
    ax = fig.add_subplot(111)

    ax.plot(displacements, angles['spice'], color='C0', label='Spice')
    # ax.axhline(y=true_angle_spice, linestyle='--', color='C0', label='Nominal Spice')
    ax.axhline(y=JF_A2_SPICE_AZI, linestyle='-.', color='C0', label='JF A2 Reco Spice')

    ax.plot(displacements, angles['rtp'], color='C1', label='RTP')
    # ax.axhline(y=true_angle_rtp, linestyle='--', color='C1', label='Nominal RTP')
    ax.axhline(y=MYL_A2_RTP_AZI, linestyle='-.', color='C1', label='MYL A2 Reco RTP')

    ax.set_xlabel('{} Displacement'.format(which_dir))

    ax.set_ylabel('Azimuth Angle')
    ax.grid()
    ax.legend()
    fig.savefig("a2_scan_{}.png".format(which_dir),edgecolor='none',bbox_inches="tight")

do_2D_scan = True
if do_2D_scan:

    xs, dx = np.linspace(-1500, 1500, 101, retstep = True)
    ys, dy = np.linspace(-1500, 1500, 101, retstep = True)

    angles_spice = np.zeros((len(xs), len(ys)))
    angles_rtp = np.zeros((len(xs), len(ys)))

    for i, x in enumerate(xs):
        for j, y in enumerate(ys):

            A2_temp = copy.copy(A2)
            A2_temp[0]+=x
            A2_temp[1]+=y
            diff_rtp = RTP - A2_temp
            diff_spice = SPIce - A2_temp
            angle_rtp = util.get_angle(diff_rtp)
            angle_spice = util.get_angle(diff_spice)
            
            angles_spice[i, j] = angle_spice
            angles_rtp[i, j] = angle_rtp

    fig, axs = plt.subplots(1, 3, figsize=(15, 5))
    
    im = axs[0].pcolormesh(xs, ys, angles_spice.T)
    cbar = plt.colorbar(im, ax=axs[0], label='Azimuthal Angle')
    axs[0].set_title('Spice Pulser')
    im.set_clim(-20, 60)
    cs = axs[0].contour(xs, ys, angles_spice.T, 
        levels=[JF_A2_SPICE_AZI], colors=('k'), linestyles=('--'))
    axs[0].clabel(cs, inline=True)


    im2 = axs[1].pcolormesh(xs, ys, angles_rtp.T)
    cbar2 = plt.colorbar(im2, ax=axs[1], label='Azimuthal Angle')
    axs[1].set_title('Rooftop Pulser ')
    im2.set_clim(-20, 60)
    cs2 = axs[1].contour(xs, ys, angles_rtp.T, 
        levels=[MYL_A2_RTP_AZI], colors=('k'), linestyles=('--'))
    axs[1].clabel(cs2, inline=True)


    cs3 = axs[2].contour(xs, ys, angles_spice.T, 
        levels=[JF_A2_SPICE_AZI], colors=('C0'), linestyles=('--'), label='Spice')
    axs[2].clabel(cs3, inline=True)
    cs4 = axs[2].contour(xs, ys, angles_rtp.T, 
        levels=[MYL_A2_RTP_AZI], colors=('C1'), linestyles=('-.'), label='RTP')
    axs[2].clabel(cs4, inline=True)
    axs[2].legend()
    
    for ax in axs:
        ax.set_xlabel('Horizontal Displacement')
        ax.set_ylabel('Vertical Displacement')
        ax.set_aspect('equal')
    plt.tight_layout()
    fig.savefig('2dscan.png')


# RTP -= A2
# SPIce -= A2
# A2-=A2
# diff = RTP
# # print("Nominal Angle {}".format(util.get_angle(diff)))


# # fig = plt.figure(figsize=(5*1.5,5*1.5))
# # ax=fig.add_subplot(111)
# # ax.plot(A2[0],A2[1], "o", label="A2")
# # ax.plot(SPIce[0],SPIce[1],"kx", label="SPIce")
# # ax.plot(RTP[0], RTP[1], 'k^', label='RTP')
# # ax.set_xlabel('Easting (km)')
# # ax.set_ylabel('Northing (km)')
# # # ax.set_xlim(-0.5,4)
# # # ax.set_ylim(-0.5,4)
# # ax.grid()
# # ax.legend(loc='lower right')
# # ax.set_aspect('equal')
# # fig.savefig("a2_spice_icecube.png",edgecolor='none',bbox_inches="tight")