# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries

tstart = "2011-08-09 07:00"
tend = "2011-08-09 10:00"
diffractionstart = parse_time('2011-08-09 07:54:39.57').datetime
diffractionend = parse_time('2011-08-09 09:03:27.57').datetime
flarestart = parse_time('2011-08-09 07:48:15.57').datetime
flareend = parse_time('2011-08-09 09:29:51.57').datetime

# tstart = "2014-09-10 17:00"
# tend = "2014-09-10 23:00"
# diffractionstart = parse_time('').datetime
# diffractionend = parse_time('').datetime
# flarestart = parse_time('2014-09-10 17:15:02.57').datetime
# flareend = parse_time('2014-09-10 18:44:38.57').datetime

tr = a.Time(tstart, tend)

results = Fido.search(tr, a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s") | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == 'SWPC'))

files = Fido.fetch(results)
goes = TimeSeries(files)

hek_results = results['hek']
flares_hek = hek_results[0]

print(type(flares_hek['fl_goescls']))
print(flares_hek['fl_goescls'])
print(type(flares_hek['fl_goescls']))
print(np.str_("Data Used"))
print(type(np.str_("Data Used")))

fig, ax = plt.subplots()
goes.plot(axes=ax)
ax.axvline(parse_time(flares_hek['event_peaktime']).datetime)
ax.axvspan(
    flarestart,
    flareend,
    alpha=0.2,
    label = np.str_("Data Used")
    # label= flares_hek['fl_goescls']
)

# SHOW RED REGION WHERE WE ARE EXCLUDING DATA
ax.axvspan(
    diffractionstart,
    diffractionend,
    alpha=0.2,
    facecolor = 'r',
    label = np.str_("Data Excluded")
)


ax.legend(loc=2)
ax.set_yscale('log')
ax.set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
ax.set_ylabel(r'SXR Flux $(\mathrm{W} \; \mathrm{m}^{-2})$')
ax.set_xlabel("Time")

plt.show()
