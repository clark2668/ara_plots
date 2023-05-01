import numpy as np
import matplotlib.pyplot as plt
import json
import pandas as pd
from scipy.interpolate import interp1d

# load the radio dec bins
gsl = open(f'data/gl_aeffs.json')
a = json.load(gsl)
radio_dec_bins = np.asarray(a['dec_edges_degrees'])

def get_bin_centers(bins):
    return (bins[1:] + bins[:-1]) * 0.5

radio_dec_bin_centers = get_bin_centers(radio_dec_bins)

def find_nearest_energy_bin(array, value):
    array = np.asarray(array)
    idx = (np.abs(array-value)).argmin()
    return idx, array[idx]

data = pd.read_csv('./data/IC86-2012-TabulatedAeff.txt', delim_whitespace=True)
data_emin = np.unique(np.asarray(data['E_min[GeV]']))
idx, en = find_nearest_energy_bin(data_emin, 1E9)
print(f"Selected energy {en:e}")

# slice to grab only the energy columns we care about
selected_data = data[data['E_min[GeV]'] == en]
dec_min = np.rad2deg(np.arccos(selected_data['cos(zenith)_min']))
dec_max = np.rad2deg(np.arccos(selected_data['cos(zenith)_min']))
dec_avg = (dec_min + dec_max)/2.
dec_avg = dec_avg-90
interpolator = interp1d(dec_avg, selected_data['Aeff[m^2]'], bounds_error=False, fill_value=0)

# to validate the interpolation
# fig = plt.figure(figsize=(7,5))
# ax = fig.add_subplot(111)
# ax.plot(
#     dec_avg,
#     selected_data['Aeff[m^2]']
# )
# radio_dec_bins = radio_dec_bins[::-1]
# weights = interpolator(radio_dec_bin_centers)
# ax.hist(
#     radio_dec_bin_centers,
#     bins=radio_dec_bins,
#     weights=weights,
#     histtype='step'
# )
# fig.tight_layout()
# fig.savefig(f'aeff_vs_zen_icecube.png')
# del fig, ax

aeff = interpolator(radio_dec_bin_centers)
weights = aeff / np.sum(aeff)

from json import JSONEncoder
class NumpyArrayEncoder(JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return JSONEncoder.default(self, obj)

radio_dec_bins = radio_dec_bins[::-1]
weights = weights[::-1]
dictionary = {
     "dec_bin_edges": np.asarray(radio_dec_bins),
     "relative_eff_areas": weights
}
json_object = json.dumps(dictionary, cls=NumpyArrayEncoder, indent=4)
with open(f"rel_areas_ic.json", "w") as outfile:
     outfile.write(json_object)