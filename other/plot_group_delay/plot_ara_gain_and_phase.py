import numpy as np
import matplotlib.pyplot as plt
from scipy import signal
import math
from scipy.fftpack import fft, ifft
from matplotlib.font_manager import FontProperties
from matplotlib import rcParams
from scipy import stats
import numpy.ma as ma

def main():
	
	#arasim=np.genfromtxt('ARA_Electronics_TotalGainPhase.txt', delimiter=',', skip_header=3, names=['FREQ','GAIN','PHASE'])
	arasim=np.genfromtxt('ARA_Electronics_TotalGain_TwoFilters.txt', delimiter=',', skip_header=3, names=['FREQ','GAIN','PHASE'])
	
	gain = arasim['GAIN']
	frequency  =  arasim['FREQ']
	phase = arasim['PHASE']
	frequency*=1e6 #convert to Hz for group delay computation
	unwrap = np.unwrap(phase)

	gdl = np.zeros_like(phase)
	max = range(frequency.shape[0])
	for i in max:
		computed_answer =0
		if (i==0 or i==256):
			computed_answer=0
		else:
			if(unwrap[i+1] ==0 or unwrap[i]==0 or unwrap[i-1]==0):
				computed_answer=0
			else:
				computed_answer = -1* (unwrap[i+1]-unwrap[i-1])/(frequency[i+1]-frequency[i-1])/(2*math.pi)
		gdl[i]=computed_answer
	gdl/=1e-9 #convert to ns for plotting
	frequency/=1e6 #convert back to MHz for plotting

	# impulse=np.array(gain,dtype='complex')
	# phases_for_fft = np.cos(phase) + 1j * np.sin(phase)
	# impulse *= phases_for_fft
	# impulse = np.fft.ifft(impulse).real

	fig = plt.figure(figsize=(18,5))
	ax1 = fig.add_subplot(1,3,1)
	ax2 = fig.add_subplot(1,3,2)
	ax3 = fig.add_subplot(1,3,3)

	labelsize=18
	
	ax1.plot(frequency,20*np.log10(arasim['GAIN']),color='steelblue',linewidth=2,label='Gain')
	ax1.set_ylabel("Voltage Gain (dB)",fontsize=labelsize)
	ax1.set_xlabel("Frequency (MHz)",fontsize=labelsize)
	ax1.set_xlim([0,1000])
	ax1.tick_params(axis='both', which='major', labelsize=16)

	# print "Max gain: ",np.max(20*np.log10(arasim['GAIN']))

	ax2.plot(frequency,unwrap,color='steelblue',linewidth=2,label='Unwrapped Phase')
	ax2.set_ylabel("Phase (rad)",fontsize=labelsize)
	ax2.set_xlabel("Frequency (MHz)",fontsize=labelsize)
	ax2.set_xlim([0,1000])
	ax2.tick_params(axis='both', which='major', labelsize=16)	

	ax3.plot(frequency,gdl,color='steelblue',linewidth=2,label='Group Delay')
	ax3.set_ylabel("Group Delay (ns)",fontsize=labelsize)
	ax3.set_xlabel("Frequency (MHz)",fontsize=labelsize)
	#ax3.set_ylim([-2,20])
	ax3.set_xlim([0,1000])
	ax3.tick_params(axis='both', which='major', labelsize=16)
	
	#ax1.grid(which='both')

	plt.tight_layout()
	fig.savefig('arasim_gain_and_phase.png',edgecolor='none',bbox_inches="tight") #save the figure

	# fig2 = plt.figure(figsize=(11.5,8.5/4))
	# ax21 = fig2.add_subplot(1,1,1)
	# time = np.arange(0,128,0.5)
	# ax21.plot(time,impulse)
	# print impulse
	# fig2.savefig("impulse_test.png")

main()
