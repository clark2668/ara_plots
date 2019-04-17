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

	fig = plt.figure(figsize=(2*5,6*5))
	num_total=6
	
	ax_a2=fig.add_subplot(num_total,1,1)
	ax_a2.hist([years_and_months],bins=bins,weights=data['A2Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	#, ,fill=False, stacked=True, histtype='step', color='red', linewidth=4, label='Testbed')
	ax_a2.set_ylabel('Fractional Uptime', size=22)

	
	ax_a3=fig.add_subplot(num_total,1,2)
	ax_a3.hist([years_and_months],bins=bins,weights=data['A3Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a3.set_ylabel('Fractional Uptime', size=22)

	ax_list = fig.axes
	for ax in ax_list:
		# ax.set_xlim([start_year-.3,stop_year+0.3]) #set the x limits of the plot
		ax.xaxis.set_ticks(np.arange(start_year, stop_year+1, 1))
		ax.set_ylim([0,1.2]) #set the x limits of the plot
		ax.set_xticklabels(['2013','2014','2015','2016','2017'], size=22)
		ax.yaxis.set_ticks([0.5,1])
		ax.set_yticklabels(['0.5','1.0'], size=22)
		ax.axhline(y=1.0, color='k', linestyle='--')
		ax.tick_params(labelsize=22)

	
	#ax5.tick_params(axis='both', which='major', labelsize=25)
	#ax5.legend(frameon=False,  fontsize=19, loc='upper left')
	fig.savefig('uptimes_a23.pdf',edgecolor='none',bbox_inches="tight")


	
main()
	
