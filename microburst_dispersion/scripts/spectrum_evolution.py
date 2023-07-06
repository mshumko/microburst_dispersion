import dateutil.parser

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

from microburst_dispersion.firebird import Hires

plot_window_s=1
time_range = (
    dateutil.parser.parse('2015-08-27T12:40:37.630000'), 
    dateutil.parser.parse('2015-08-27T12:40:37.800000')
)
max_channel = 4

fb_id = 3

hr = Hires(fb_id, time_range[0]).load()
idx = np.where((hr['Time'] >= time_range[0]) & (hr['Time'] <= time_range[1]))[0]
time = hr['Time'][idx]
col_flux = hr['Col_flux'][idx, :max_channel]
center_energy = np.array([
    float(s.split()[0].replace('>', '')) 
    for s in np.array(hr.attrs['Col_counts']['ELEMENT_LABELS'][:max_channel])
    ])


fig, ax = plt.subplots(figsize=(5, 5))
cmap = matplotlib.colormaps['viridis']
colors = ['k', 'purple', 'b', 'c', 'g', 'r']
# cmap = cmap.reversed()
for i, (time_i, col_flux_i) in enumerate(zip(time[::2], col_flux[::2, :])):
    ax.plot(center_energy, col_flux_i, label=f'{time_i:%S.%f}', color=colors[i])

ax.legend(title='Time [Seconds]')
plt.suptitle(f'FU{fb_id} Energy Spectrum\n{str(time_range[0])[:-3]} to {str(time_range[1])[:-3]}')
ax.set(yscale='log', xlabel='Energy [keV]', ylabel='Flux')
plt.tight_layout()
plt.show()
pass