#!/bin/python

import sys
import matplotlib.pyplot as plt
import os
import numpy as np

# def read_efficiency_file(STATION, TYPE):
#    data = []
#    '''   
#    if TYPE in ['2','4']:
#       fp = STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_tunedPostCutCoherenceThermalCut_snrCut_signalEfficiency.csv'
#    else:
#       fp = STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_coherenceThermalCut_snrCut_signalEfficiency.csv'
#    '''
#    fp=STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_tunedPostCutThermalSNRSurfaceCut_signalEfficiency.csv'
#    #fp = STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_coherenceThermalCut_snrCut_signalEfficiency.csv'
#    #with open(STATION+'_type'+TYPE+'_signalEfficiency_fullCWCut_thermalImpCut.csv') as f:
#    #with open(STATION+'_type'+TYPE+'_impulsivityFilter1_signalEfficiency.csv') as f:
#    #with open(STATION+'_type'+TYPE+'_impulsivityFilter1_noMaskSat_flattenSat_signalEfficiency.csv') as f:
#    #with open(STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_newCWCut_signalEfficiency.csv') as f:
#    #with open(STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_coherenceThermalCut_snrCut_signalEfficiency.csv') as f:
#    with open(fp) as f:
#       for row in f:
#          data.append(map(float, row.rstrip().split(',')))
#    #map(float, data)      
#    return data

# def get_script_path():
#    return os.path.dirname(os.path.realpath(sys.argv[0]))


# def sum_tables(t1, t2):
#    nRow = len(t1)
#    nCol = len(t1[0])
#    for i in range(nRow):
#       for j in range(1,nCol):
#          t1[i][j] += t2[i][j]
#    return t1


# STATION = sys.argv[1]
# #TYPE = sys.argv[2]
# eff_types = []
# types = ['1','2','3','4','5']

# for TYPE in types:
#    eff_types.append( read_efficiency_file(STATION, TYPE) )

# eff = sum_tables(eff_types[0], sum_tables(eff_types[1], sum_tables(eff_types[2], sum_tables(eff_types[3], eff_types[4]))))

# '''
# len_rows = len(eff_types[0])
# len_cols = len(eff_types[0][0])
# print ('len_rows: %d') % (len_rows)

# eff = []

# for TYPE in range(len(types)):
#    for row in range(len_rows):
#       for col in range(len_cols):
#          if row == 0 and col == 0:
# '''



# energy = []
# trig = []
# offset = []
# nchnl = []
# #imp = []
# #sat = []
# cw = []
# #cw_new = []
# #dp = []
# thermal = []
# #thermal_imp = []
# snr = []
# cal = []
# surface = []
# #noisy = []

# row_count = 0
# for row in eff:
#    energy.append(float(row[0])-9)
#    trig.append(float(row[1])/float(row[1]))
#    offset.append(float(row[2])/float(row[1]))
#    nchnl.append(float(row[3])/float(row[1]))
#    #imp.append(float(row[3])/float(row[1]))
#    #sat.append(float(row[4])/float(row[1]))
#    cw.append(float(row[4])/float(row[1]))
#    #cw_new.append(float(row[5])/float(row[1]))
#    #dp.append(float(row[5])/float(row[1]))
#    thermal.append(float(row[5])/float(row[1]))
#    #thermal_imp.append(float(row[6])/float(row[1]))
#    snr.append(float(row[6])/float(row[1]))
#    cal.append(float(row[7])/float(row[1]))
#    surface.append(float(row[8])/float(row[1]))
#    #noisy.append(float(row[9])/float(row[1]))
#    row_count+=1


plt.figure(figsize=(2*11,2*8.5))
ax = plt.subplot(2,2,1)
ax2 = plt.subplot(2,2,2)
'''
ax.bar(energy, trig, 0.5, align='center', label = 'trig', color='black')
ax.bar(energy, offset, 0.5, align='center', label='offset block', color='blue')
ax.bar(energy, nchnl, 0.5, align='center', label='n-chan', color='green')
ax.bar(energy, cw_new, 0.5, align='center', label='CW', color='red')
#ax.bar(energy, dp, align='center', label='deep pulser', color='magenta')
ax.bar(energy, thermal, 0.5, align='center', label='thermal', color='yellow')
ax.bar(energy, cal, 0.5, align='center', label='calpulser', color='cyan')
ax.bar(energy, surface, 0.5, align='center', label='surface', color='gray')
'''

# ax.plot(energy, trig, '.k-', label = 'Trigger', lw=5)
# ax.plot(energy, offset, '.b-', label='Offset Block', lw=5)
# ax.plot(energy, nchnl, 'black', label='Event Filter', lw=5) #n-chan
#ax.plot(energy, imp, '.b-', label='Impulsivity Filter (0.001)', lw=5)
#ax.plot(energy, sat, color='pink', linestyle='-', marker='.', label='saturation', )
# ax.plot(energy, cw, '.r--', label='Low Freq. Cut', lw=5)
#ax.plot(energy, cw_new, '.r-.', label='CW')
#ax.bar(energy, dp, align='center', label='deep pulser', color='magenta')
# ax.plot(energy, thermal, linestyle='-', marker='.', color='cyan', label='Thermal Cut', lw=5)
#ax.plot(energy, thermal_imp, color='orange', linestyle='-.', marker='.', label='thermal imp.', lw=5)
# ax.plot(energy, snr, color='red', linestyle='--', marker='.', label='+2D Cut', lw=5) #SNR cut
# ax.plot(energy, cal, linestyle='-', marker='.', color='green', label='+Calpulser Cut', lw=5)
# ax.plot(energy[2:-1], surface[2:-1], color='red', linestyle='--', label='Analysis A', lw=5)
# print surface
#ax.bar(energy, noisy, align='center', label='noisy run', color='orange')

# and now for Brian's part
eff_ebins = np.array([17.25, 17.75, 18.25, 18.75, 19.25, 19.75, 20.25])
eff_ebins-=9 #convert to GeV
eff_c1 = np.array([0.34, 0.42, 0.48, 0.56, 0.58, 0.66, 0.71])
eff_c2 = np.array([0.25, 0.34, 0.40, 0.44, 0.48, 0.60, 0.64])
eff_c3 = np.array([0.33, 0.43, 0.49, 0.58, 0.63, 0.67, 0.87])
eff_c4 = np.array([0.29, 0.38, 0.46, 0.55, 0.58, 0.65, 0.69])
eff_c5 = np.array([0.31, 0.40, 0.46, 0.54, 0.60, 0.67, 0.80])
frac_uptime = np.array([0.165, 0.132, 0.087, 0.375, 0.242])
eff_ave = (frac_uptime[0]*eff_c1) + (frac_uptime[1]*eff_c2) + (frac_uptime[2]*eff_c3) + (frac_uptime[3]*eff_c4) + (frac_uptime[4]*eff_c5)
ax.plot(eff_ebins, eff_ave, color='blue', linestyle='-.', label='This Work', linewidth=7)

data_tb = np.genfromtxt("testbed_efficiency.csv",delimiter=',',skip_header=1,names=['energy_logev','eff'])
tb_logeV = data_tb['energy_logev']
tb_logeV-=9
tb_eff = data_tb['eff']
ax.plot(tb_logeV[1:-6], tb_eff[1:-6], '-', color='grey', label='Testbed Analysis', linewidth=5)


brian_data = np.genfromtxt("brian_eff_vs_snr.csv",delimiter=',',skip_header=1,names=['snr','c1','c2','c3','c4','c5'])
snr_bins = brian_data['snr']
eff_vs_snr_c1 = brian_data['c1']
eff_vs_snr_c2 = brian_data['c2']
eff_vs_snr_c3 = brian_data['c3']
eff_vs_snr_c4 = brian_data['c4']
eff_vs_snr_c5 = brian_data['c5']

brian_eff_vs_snr_ave = (frac_uptime[0]*eff_vs_snr_c1) + (frac_uptime[1]*eff_vs_snr_c2) + (frac_uptime[2]*eff_vs_snr_c3) + (frac_uptime[3]*eff_vs_snr_c4) + (frac_uptime[4]*eff_vs_snr_c5)
ax2.plot(snr_bins,brian_eff_vs_snr_ave, label='This Work',linestyle='-.',linewidth=7, color='blue')


# myl_data = np.genfromtxt("myl_eff_vs_snr.csv",delimiter=',',skip_header=1,names=['snr','c1','c2','c3','c4','c5'])
# snr_bins = myl_data['snr']
# eff_vs_snr_c1 = myl_data['c1']
# eff_vs_snr_c2 = myl_data['c2']
# eff_vs_snr_c3 = myl_data['c3']
# eff_vs_snr_c4 = myl_data['c4']
# eff_vs_snr_c5 = myl_data['c5']
# myl_eff_vs_snr_ave = (frac_uptime[0]*eff_vs_snr_c1) + (frac_uptime[1]*eff_vs_snr_c2) + (frac_uptime[2]*eff_vs_snr_c3) + (frac_uptime[3]*eff_vs_snr_c4) + (frac_uptime[4]*eff_vs_snr_c5)
# ax2.plot(snr_bins,myl_eff_vs_snr_ave, label='Analysis A',linestyle='--',linewidth=5, color='red')

tb_data = np.genfromtxt("tb_eff_vs_snr.csv",delimiter=',',skip_header=1,names=['snr','eff'])
snr_bins = tb_data['snr']
tb_eff = tb_data['eff']
ax2.plot(snr_bins,tb_eff, label='Testbed Analysis',linestyle='-',linewidth=5, color='grey')





'''
x_label = ['trig', 'offset', 'nchnl', 'CW', 'thermal', 'thermal_imp', 'cal', 'surface']
pass_rate = []
colors = ['gray', 'blue', 'magenta', 'green', 'red', 'gold', 'cyan', 'black', 'orange', 'pink', 'teal']
for row in range(row_count):
  pass_rate.append([trig[row], offset[row], nchnl[row], cw_new[row], thermal[row], cal[row], surface[row]])
  ax.plot(range(1,8), pass_rate[row], color=colors[row], linestyle='-', marker='.', label='log(E/eV)='+str(16+row*0.5))
'''

#ax.semilogy()

#ax.set_yscale('log')
ax.set_xlim(17-9,20.5-9)
ax.set_ylim(2e-2,1.1)
ax.set_xlabel('log(E/GeV)', fontsize=30)
ax.set_ylabel('Signal Efficiency', fontsize=30)
ax.set_title('A2 Efficiency vs Energy', fontsize=30)
ax.grid(which='both',axis='y')
ax.yaxis.set_tick_params(labelsize=26)
ax.xaxis.set_tick_params(labelsize=26)


ax2.set_ylim(2e-2,1.1)
ax2.set_xlabel('SNR', fontsize=30)
ax2.set_ylabel('Signal Efficiency', fontsize=30)
ax2.set_title('A2 Efficiency vs SNR', fontsize=30)
ax2.grid(which='both',axis='y')
ax2.yaxis.set_tick_params(labelsize=26)
ax2.xaxis.set_tick_params(labelsize=26)


'''
ax.set_yscale('log')
ax.set_xlim(0,8)
ax.set_ylim(3e-2, 2)
ax.set_xticks(range(1,8))
#ax.set_yticklabels(fontsize=14)
ax.yaxis.set_tick_params(labelsize=14)
ax.set_xticklabels(x_label, fontsize=18)
ax.set_ylabel('signal efficiency', fontsize=18)
ax.set_title(STATION+' config '+TYPE+' signal efficieny')
ax.grid(which='both',axis='y')
'''
ax.legend(loc='upper left', fontsize=26)
# ax2.legend(loc='upper left', fontsize=26)
plt.tight_layout()
# PWD = get_script_path()
#plt.savefig(PWD+'/'+STATION+'_type'+TYPE+'_signalEfficiencyLineAlongCuts.png', bbox_inches='tight')
#plt.savefig(PWD+'/'+STATION+'_type'+TYPE+'_signalEfficiency_impulsivityFilter1.png', bbox_inches='tight')
#plt.savefig(PWD+'/'+STATION+'_signalEfficiency_impulsivityFilter1_noMaskSat_flattenSat.png', bbox_inches='tight')
#plt.savefig(PWD+'/'+STATION+'_signalEfficiency_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_tunedPostCutThermalSNRSurfaceCut.png', bbox_inches='tight')
plt.savefig("A2_efficiency.png",dpi=300,edgecolor='none',bbox_inches="tight")
# plt.show()