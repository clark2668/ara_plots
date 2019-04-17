#!/bin/python

import sys
import matplotlib.pyplot as plt
import os

def read_efficiency_file(STATION, TYPE):
   data = []
   '''   
   if TYPE in ['2','4']:
      fp = STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_tunedPostCutCoherenceThermalCut_snrCut_signalEfficiency.csv'
   else:
      fp = STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_coherenceThermalCut_snrCut_signalEfficiency.csv'
   '''
   fp=STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_tunedPostCutThermalSNRSurfaceCut_signalEfficiency.csv'
   #fp = STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_coherenceThermalCut_snrCut_signalEfficiency.csv'
   #with open(STATION+'_type'+TYPE+'_signalEfficiency_fullCWCut_thermalImpCut.csv') as f:
   #with open(STATION+'_type'+TYPE+'_impulsivityFilter1_signalEfficiency.csv') as f:
   #with open(STATION+'_type'+TYPE+'_impulsivityFilter1_noMaskSat_flattenSat_signalEfficiency.csv') as f:
   #with open(STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_newCWCut_signalEfficiency.csv') as f:
   #with open(STATION+'_type'+TYPE+'_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_coherenceThermalCut_snrCut_signalEfficiency.csv') as f:
   with open(fp) as f:
      for row in f:
         data.append(map(float, row.rstrip().split(',')))
   #map(float, data)      
   return data

def get_script_path():
   return os.path.dirname(os.path.realpath(sys.argv[0]))


def sum_tables(t1, t2):
   nRow = len(t1)
   nCol = len(t1[0])
   for i in range(nRow):
      for j in range(1,nCol):
         t1[i][j] += t2[i][j]
   return t1


STATION = sys.argv[1]
#TYPE = sys.argv[2]
eff_types = []
types = ['1','2','3','4','5']

for TYPE in types:
   eff_types.append( read_efficiency_file(STATION, TYPE) )

eff = sum_tables(eff_types[0], sum_tables(eff_types[1], sum_tables(eff_types[2], sum_tables(eff_types[3], eff_types[4]))))

'''
len_rows = len(eff_types[0])
len_cols = len(eff_types[0][0])
print ('len_rows: %d') % (len_rows)

eff = []

for TYPE in range(len(types)):
   for row in range(len_rows):
      for col in range(len_cols):
         if row == 0 and col == 0:
'''



energy = []
trig = []
offset = []
nchnl = []
#imp = []
#sat = []
cw = []
#cw_new = []
#dp = []
thermal = []
#thermal_imp = []
snr = []
cal = []
surface = []
#noisy = []

row_count = 0
for row in eff:
   energy.append(float(row[0])-9)
   trig.append(float(row[1])/float(row[1]))
   offset.append(float(row[2])/float(row[1]))
   nchnl.append(float(row[3])/float(row[1]))
   #imp.append(float(row[3])/float(row[1]))
   #sat.append(float(row[4])/float(row[1]))
   cw.append(float(row[4])/float(row[1]))
   #cw_new.append(float(row[5])/float(row[1]))
   #dp.append(float(row[5])/float(row[1]))
   thermal.append(float(row[5])/float(row[1]))
   #thermal_imp.append(float(row[6])/float(row[1]))
   snr.append(float(row[6])/float(row[1]))
   cal.append(float(row[7])/float(row[1]))
   surface.append(float(row[8])/float(row[1]))
   #noisy.append(float(row[9])/float(row[1]))
   row_count+=1


plt.figure(figsize=(12,12))
ax = plt.subplot(1,1,1)
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
ax.plot(energy, nchnl, 'black', label='Event Filter', lw=5) #n-chan
#ax.plot(energy, imp, '.b-', label='Impulsivity Filter (0.001)', lw=5)
#ax.plot(energy, sat, color='pink', linestyle='-', marker='.', label='saturation', )
# ax.plot(energy, cw, '.r--', label='Low Freq. Cut', lw=5)
#ax.plot(energy, cw_new, '.r-.', label='CW')
#ax.bar(energy, dp, align='center', label='deep pulser', color='magenta')
# ax.plot(energy, thermal, linestyle='-', marker='.', color='cyan', label='Thermal Cut', lw=5)
#ax.plot(energy, thermal_imp, color='orange', linestyle='-.', marker='.', label='thermal imp.', lw=5)
ax.plot(energy, snr, color='red', linestyle='--', marker='.', label='+2D Cut', lw=5) #SNR cut
ax.plot(energy, cal, linestyle='-', marker='.', color='green', label='+Calpulser Cut', lw=5)
ax.plot(energy, surface, color='blue', linestyle='--', marker='.', label='+Surface Cut', lw=5)
print surface
#ax.bar(energy, noisy, align='center', label='noisy run', color='orange')

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
ax.set_xlim(15-9,22-9)
ax.set_ylim(2e-2,1.1)
ax.set_xlabel('log(E/eV)', fontsize=30)
ax.set_ylabel('Signal Efficiency', fontsize=30)
#ax.set_title(STATION+' Config '+TYPE+' Signal Efficieny')
ax.set_title(STATION+' Signal Efficiency', fontsize=30)
ax.grid(which='both',axis='y')
ax.yaxis.set_tick_params(labelsize=26)
ax.xaxis.set_tick_params(labelsize=26)

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
plt.legend(loc='lower right', fontsize=26)
plt.tight_layout()
PWD = get_script_path()
#plt.savefig(PWD+'/'+STATION+'_type'+TYPE+'_signalEfficiencyLineAlongCuts.png', bbox_inches='tight')
#plt.savefig(PWD+'/'+STATION+'_type'+TYPE+'_signalEfficiency_impulsivityFilter1.png', bbox_inches='tight')
#plt.savefig(PWD+'/'+STATION+'_signalEfficiency_impulsivityFilter1_noMaskSat_flattenSat.png', bbox_inches='tight')
#plt.savefig(PWD+'/'+STATION+'_signalEfficiency_vnchnl3NoMasking_noMaskSat_flattenSat_snrMode1_tunedPostCutThermalSNRSurfaceCut.png', bbox_inches='tight')
plt.savefig("A2_efficiency.png")
plt.show()

