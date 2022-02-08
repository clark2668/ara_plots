import sys
import matplotlib.pyplot as plt
import os
import numpy as np

fig = plt.figure()
ax = plt.subplot(111)

data = np.genfromtxt('baseline_fractions.csv',
    delimiter=',',
    skip_header=1,
    names=['energy', 
        'na1', 'na2', 'na3', 'na4', 'na5', 'na6', 'na7', 'na8',
        'frac_hybrid', 'frac_shallow', 'frac_coinc'
        ]
)

ax.plot(data['energy'], data['frac_hybrid'], #'o-', 
    linestyle='-', color='C2',
    linewidth=4, label='Hybrid Stations')
ax.plot(data['energy'], data['frac_shallow'], #'s--', 
    linestyle='--', color='C1',
    linewidth=4, label='Surface Stations')
ax.plot(data['energy'], data['frac_coinc'], #'^:', 
    linestyle=':', color='C0',
    linewidth=4, label='Coincidences ("Golden Events")')

ax.set_xlabel(r'$log_{10}(E_{\nu})$ (eV)', fontsize=20)
ax.set_ylabel(r'Fraction  of  $V_{eff}$', fontsize=20)
ax.set_ylim([0, 1])
ax.yaxis.set_tick_params(labelsize=20)
ax.xaxis.set_tick_params(labelsize=20)
plt.tight_layout()
fig.savefig('coincidence_fractions.png')
