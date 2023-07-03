"""
Calculate how large the ducted region need to be if the TOF dispersion is due
to the difference in how chorus waves propagate within the duct.
"""
import dateutil.parser

import IRBEM
import numpy as np

import microburst_dispersion.firebird

Re = 6378.14  # km

time_range = (
    dateutil.parser.parse('2015-08-27T12:40:37.5'), 
    dateutil.parser.parse('2015-08-27T12:40:37.7')
    )
fb_id = 3

hr = microburst_dispersion.firebird.Hires(fb_id, time_range[0]).load()
xgeo = np.full((2, 3), np.nan)

model = IRBEM.MagFields()
for i, time in enumerate(time_range):
    hr_idx = np.where(hr['Time'] > time)[0][0]
    assert np.abs(hr['Time'][hr_idx] - time).total_seconds() < hr.attrs['CADENCE']
    X = {'time':time, 'x1':hr['Alt'][hr_idx], 'x2':hr['Lat'][hr_idx], 'x3':hr['Lon'][hr_idx]}
    xgeo[i, :] = model.find_magequator(X, {'Kp':20})['XGEO']

diff = xgeo[0, :] - xgeo[1, :]
distance = Re*np.linalg.norm(diff)
print(f"The microburst observation during FU{fb_id}'s orbit corresponded to an "
      f"{round(distance)} km equatorial size (mostly radial).")