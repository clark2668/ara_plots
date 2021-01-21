import numpy as np
import matplotlib.pyplot as plt
import itertools


def main():
	
	data = np.genfromtxt('livetime_frac_at12345.csv',delimiter=',',skip_header=1,
		names=['Year','Month','Date','TestbedFrac','A1Frac','A2Frac','A3Frac','A4Frac','A5Frac'])
	years = data['Year']
	months = data['Month']
	years_and_months = years + (months/12.) - 1./12.
	bins = years_and_months
	print(bins)

	start_year=2012
	stop_year=2021
	# bins = 49

	fig = plt.figure(figsize=(5,10))
	num_total=5
	
	horz_shift=0.225

	ax_a1=fig.add_subplot(num_total,1,1)
	ax_a1.hist([years_and_months],bins=bins,weights=data['A1Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a1.set_ylabel('A1 Uptime (%)', size=11)
	ax_a1.text(horz_shift, 0.90, '1598 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a1.transAxes,size=11)

	ax_a2=fig.add_subplot(num_total,1,2)
	ax_a2.hist([years_and_months],bins=bins,weights=data['A2Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a2.set_ylabel('A2 Uptime (%)', size=11)
	ax_a2.text(horz_shift, 0.90, '2219 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a2.transAxes,size=11)

	ax_a3=fig.add_subplot(num_total,1,3)
	ax_a3.hist([years_and_months],bins=bins,weights=data['A3Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a3.set_ylabel('A3 Uptime (%)', size=11)
	ax_a3.text(horz_shift, 0.90, '2103 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a3.transAxes,size=11)

	ax_a4=fig.add_subplot(num_total,1,4)
	ax_a4.hist([years_and_months],bins=bins,weights=data['A4Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a4.set_ylabel('A4 Uptime (%)', size=11)
	ax_a4.text(horz_shift, 0.90, '1031 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a4.transAxes,size=11)

	ax_a5=fig.add_subplot(num_total,1,5)
	ax_a5.hist([years_and_months],bins=bins,weights=data['A5Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
	ax_a5.set_ylabel('A5 Uptime (%)', size=11)
	ax_a5.text(horz_shift, 0.90, '729 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a5.transAxes,size=11)

	ax_list = fig.axes
	for ax in ax_list:
		ax.set_xlim([start_year,2021]) #set the x limits of the plot
		ax.xaxis.set_ticks(np.array([2012,2013,2014,2015,2016,2017,2018,2019,2020,2021]))
		# ax.xaxis.set_ticks(np.arange(start_year, stop_year+1, 1))
		ax.set_xticklabels(['','','','','','','','',''])
		
		ax.set_ylim([0,1.2]) #set the y limits of the plot
		ax.yaxis.set_ticks([0.5,1])
		ax.set_yticklabels(['0.5','1.0'], size=12)
		
		ax.axhline(y=1.0, color='k', linestyle='--')
		ax.tick_params(labelsize=12)
		ax.grid(True,'major','x',linestyle='--');


	# ax_a1.get_xaxis().set_visible(False)
	ax_a5.set_xlabel('Year',size=12)
	# ax_a5.xaxis.set_ticks(np.arange(start_year, stop_year+1, 1))
	# ax_a5.xaxis.set_ticks(np.arange(2012, 2020+1, 1))
	ax_a5.set_xticklabels(['2012','\'13','\'14','\'15','\'16','\'17','\'18','\'19','\'20','\'21'],size=11)
	
	#ax5.tick_params(axis='both', which='major', labelsize=25)
	#ax5.legend(frameon=False,  fontsize=19, loc='upper left')
	fig.savefig('livetime_a12345.png',edgecolor='none',bbox_inches="tight",dpi=300, w_pad=1.5)

main()


# def main_2column():
	
# 	data = np.genfromtxt('livetime_frac_at12345.csv',delimiter=',',skip_header=1,
# 		names=['Year','Month','Date','TestbedFrac','A1Frac','A2Frac','A3Frac','A4Frac','A5Frac'])
# 	years = data['Year']
# 	months = data['Month']
# 	years_and_months = years + (months/12.) - 1./12.
# 	bins = years_and_months

# 	start_year=2012
# 	stop_year=2021
# 	# bins = 49

# 	fig = plt.figure(figsize=(5,10))
# 	num_total=5
	
# 	horz_shift=0.225

# 	fig, ((ax_a1, ax_a2,ax_a3), (ax_a4, ax_a5, ax_a6)) = plt.subplots(2, 3)

# 	# ax_a1=fig.add_subplot(num_total,1,1)
# 	ax_a1.hist([years_and_months],bins=bins,weights=data['A1Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
# 	ax_a1.set_ylabel('A1 Uptime (%)', size=11)
# 	ax_a1.text(horz_shift, 0.90, '1353 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a1.transAxes,size=11)

# 	# ax_a2=fig.add_subplot(num_total,1,2)
# 	ax_a2.hist([years_and_months],bins=bins,weights=data['A2Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
# 	ax_a2.set_ylabel('A2 Uptime (%)', size=11)
# 	ax_a2.text(horz_shift, 0.90, '1989 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a2.transAxes,size=11)

# 	# ax_a3=fig.add_subplot(num_total,1,3)
# 	ax_a3.hist([years_and_months],bins=bins,weights=data['A3Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
# 	ax_a3.set_ylabel('A3 Uptime (%)', size=11)
# 	ax_a3.text(horz_shift, 0.90, '1858 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a3.transAxes,size=11)

# 	# ax_a4=fig.add_subplot(num_total,2,1)
# 	ax_a4.hist([years_and_months],bins=bins,weights=data['A4Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
# 	ax_a4.set_ylabel('A4 Uptime (%)', size=11)
# 	ax_a4.text(horz_shift, 0.90, '786 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a4.transAxes,size=11)

# 	# ax_a5=fig.add_subplot(num_total,2,2)
# 	ax_a5.hist([years_and_months],bins=bins,weights=data['A5Frac'],range=(start_year,stop_year),edgecolor='black',color='darkseagreen',histtype='stepfilled')
# 	ax_a5.set_ylabel('A5 Uptime (%)', size=11)
# 	ax_a5.text(horz_shift, 0.90, '538 Days', horizontalalignment='right', verticalalignment='center', transform=ax_a5.transAxes,size=11)

# 	ax_list = fig.axes
# 	for ax in ax_list:
# 		ax.set_xlim([start_year,stop_year]) #set the x limits of the plot
# 		ax.xaxis.set_ticks(np.arange(start_year, stop_year+1, 1))
# 		ax.set_xticklabels(['','','','','','','','',''])
		
# 		ax.set_ylim([0,1.2]) #set the y limits of the plot
# 		ax.yaxis.set_ticks([0.5,1])
# 		ax.set_yticklabels(['0.5','1.0'], size=12)
		
# 		ax.axhline(y=1.0, color='k', linestyle='--')
# 		ax.tick_params(labelsize=12)
# 		ax.grid(True,'major','x',linestyle='--');


# 	# ax_a1.get_xaxis().set_visible(False)
# 	ax_a5.set_xlabel('Year',size=12)
# 	ax_a5.xaxis.set_ticks(np.arange(start_year, stop_year+1, 1))
# 	ax_a5.set_xticklabels(['2012','2013','2014','2015','2016','2017','2018','2019','2020','2021'],size=11)
	
# 	#ax5.tick_params(axis='both', which='major', labelsize=25)
# 	#ax5.legend(frameon=False,  fontsize=19, loc='upper left')
# 	fig.savefig('livetime_a12345.png',edgecolor='none',bbox_inches="tight",dpi=300, w_pad=1.5)

# main_2column()
	

	
