import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('livetime_vs_time.csv',delimiter=',',skip_header=1,
	names=['Year','Month','Date','Sum'])
date = data['Date']
livetime = data['Sum']
start_year=2012
stop_year=2021

fig = plt.figure(figsize=(7,5 ))
ax=fig.add_subplot(1,1,1)
ax.plot(date,livetime,'o')
ax.tick_params(labelsize=14)
ax.grid(True,'major','x',linestyle='--');
ax.set_xlabel('Year',size=14)
ax.set_ylabel(r'Total Accumulated Livetime [station-years]', size=14)
plt.tight_layout()
fig.savefig('livetime_vs_time.png',edgecolor='none',bbox_inches="tight",dpi=300, w_pad=1.5)
