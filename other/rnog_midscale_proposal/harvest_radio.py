import numpy as np
import matplotlib.pyplot as plt
import json


def find_nearest_energy_bin(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

def get_bin_centers(bins):
    return (bins[1:] + bins[:-1]) * 0.5

def get_bin_widths(bins):
     return bins[:-1] - bins[1:]

which = 'sp'

gsl = open(f'data/{which}_aeffs.json')
a = json.load(gsl)

# dict_keys(['dec_edges_degrees', 'effective_areas_m^2_bands', 'energies_eV'])

energies = a['energies_eV']
# idx, en = find_nearest_energy_bin(energies, 1.05e18)
idx, en = find_nearest_energy_bin(energies, 31E15)
print(f"Picked energy: {en}")

aeff = []
# loop over all declination bins and pull out
# the effective area at the energy we care about
for key in a['effective_areas_m^2_bands'].keys():
        vals = a['effective_areas_m^2_bands'][key]
        this_eff_at_e = vals[idx]
        aeff.append(this_eff_at_e)

aeff = np.asarray(aeff)

dec_bin_edges = np.asarray(a['dec_edges_degrees'])
sindec_bin_edges = np.sin(np.deg2rad(dec_bin_edges))
dsintheta = get_bin_widths(sindec_bin_edges)
weights = aeff/np.sum(aeff)

dec_bin_edges = dec_bin_edges[::-1] # reverse order
weights = weights[::-1]
print(f"Weights {np.sum(weights)}")

dec_bin_centers = get_bin_centers(dec_bin_edges) # for plotting

fig = plt.figure(figsize=(7,5))
ax = fig.add_subplot(111)

ax.hist(
    dec_bin_centers,
    bins=dec_bin_edges,
    weights=weights,
    histtype='step',
)

ax.set_ylim([0, 0.3])
ax.set_xlim([-90, 90])
# ax.set_yscale('log')
ax.grid()
fig.tight_layout()
fig.savefig(f'aeff_vs_zen_{idx}_{which}.png')
del fig, ax

from json import JSONEncoder
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

dictionary = {
     "dec_bin_edges": np.asarray(dec_bin_edges),
     "relative_eff_areas": weights
}
json_object = json.dumps(dictionary, cls=NumpyArrayEncoder, indent=4)
with open(f"rel_areas_{which}.json", "w") as outfile:
     outfile.write(json_object)
