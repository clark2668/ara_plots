import numpy as np
import matplotlib.pyplot as plt
import itertools


def main():
	
	data = np.genfromtxt('a23_livetime.csv',delimiter=',',skip_header=1,names=['Year','Month','A2Frac','A3Frac'])
	years = data['Year']
	months = data['Month']
	years_and_months = years + (months/12.) - 1./12.
	bins = years_and_months

	start_year=2013
	stop_year=2017
	# bins = 49

	fig = plt.figure(figsize=(10,5))
	num_total=2
	
	ax_a2=fig.add_subplot(num_total,2,1)
	ax_a2.hist([years_and_months],bins=bins,weights=data['A2Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	#, ,fill=False, stacked=True, histtype='step', color='red', linewidth=4, label='Testbed')
	ax_a2.set_ylabel('A2 Fractional Livetime', size=11)
	ax_a2.set_xlabel('Year', size=11)
	ax_a2.text(0.99, 0.90, '1142 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a2.transAxes,size=11)

	
	ax_a3=fig.add_subplot(num_total,2,2)
	ax_a3.hist([years_and_months],bins=bins,weights=data['A3Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a3.set_ylabel('A3 Fractional Livetime', size=11)
	ax_a3.set_xlabel('Year', size=11)
	ax_a3.text(0.99, 0.90, '1078 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a3.transAxes,size=11)


	ax_list = fig.axes
	for ax in ax_list:
		# ax.set_xlim([start_year-.3,stop_year+0.3]) #set the x limits of the plot
		ax.xaxis.set_ticks(np.arange(start_year, stop_year+1, 1))
		ax.set_ylim([0,1.2]) #set the x limits of the plot
		ax.set_xticklabels(['2013','2014','2015','2016','2017'], size=12)
		ax.yaxis.set_ticks([0.5,1])
		ax.set_yticklabels(['0.5','1.0'], size=12)
		ax.axhline(y=1.0, color='k', linestyle='--')
		ax.tick_params(labelsize=12)

	
	#ax5.tick_params(axis='both', which='major', labelsize=25)
	#ax5.legend(frameon=False,  fontsize=19, loc='upper left')
	fig.savefig('livetimes_a23.png',edgecolor='none',bbox_inches="tight",dpi=300, w_pad=1.5)


	
main()
	
