import sys
import matplotlib.pyplot as plt
import os
import numpy as np

#first eff vs SNR
fig = plt.figure(figsize=(10,10))
ax = plt.subplot(1,1,1)

data_a2_snr = np.genfromtxt("A2_signalEfficiencyVsChannelInWindowSNR_E-2.13.csv",delimiter=',',skip_header=0,names=['snr','eff'])
a2_snr = data_a2_snr['snr']
a2_snr_eff = data_a2_snr['eff']
ax.plot(a2_snr, a2_snr_eff, '-', color='blue', label='A2, 4 years', linewidth=5)

data_tb_snr = np.genfromtxt("testbed_eff_vs_snr.csv",delimiter=',',skip_header=1,names=['snr','eff'])
tb_snr = data_tb_snr['snr']
tb_snr_eff = data_tb_snr['eff']
ax.plot(tb_snr, tb_snr_eff, '--', color='grey', label='Prototype', linewidth=5)


ax.set_ylim(0,1.1)
ax.set_xlim(0,20)
ax.set_xlabel('SNR', fontsize=30)
ax.set_ylabel('Signal Efficiency', fontsize=30)
# ax.set_title('A2 Efficiency vs SNR', fontsize=30)
ax.grid(which='both',axis='y')
ax.yaxis.set_tick_params(labelsize=26)
ax.xaxis.set_tick_params(labelsize=26)
ax.legend(loc='upper left', fontsize=26)

plt.tight_layout()
plt.savefig("eff_vs_snr_comparison.png",dpi=300,edgecolor='none',bbox_inches="tight")



#first eff vs energy
fig2 = plt.figure(figsize=(10,10))
ax2 = plt.subplot(1,1,1)

data_a2_energy = np.genfromtxt("ARA02_liveTimeWeightedSignalEfficiency.csv",delimiter=',',skip_header=0,names=['energy','eff','err'])
a2_energy = data_a2_energy['energy']
a2_energy_eff = data_a2_energy['eff']
a2_energy = np.log10(a2_energy)
a2_energy+=9
ax2.plot(a2_energy, a2_energy_eff, '-', color='blue', label='A2, 4 years', linewidth=5)

data_tb_energy = np.genfromtxt("testbed_eff_vs_energy.csv",delimiter=',',skip_header=1,names=['energy','eff'])
tb_energy = data_tb_energy['energy']
tb_energy_eff = data_tb_energy['eff']
ax2.plot(tb_energy, tb_energy_eff, '--', color='grey', label='Prototype', linewidth=5)


ax2.set_ylim(0,1.1)
ax2.set_xlim(15.5,21.5)
ax2.set_xlabel(r'Energy [log10(eV)]', fontsize=30)
ax2.set_ylabel('Signal Efficiency', fontsize=30)
# ax.set_title('A2 Efficiency vs SNR', fontsize=30)
ax2.grid(which='both',axis='y')
ax2.yaxis.set_tick_params(labelsize=26)
ax2.xaxis.set_tick_params(labelsize=26)
ax2.legend(loc='upper left', fontsize=26)

plt.tight_layout()
plt.savefig("eff_vs_energy_comparison.png",dpi=300,edgecolor='none',bbox_inches="tight")