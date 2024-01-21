# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

import matplotlib.pyplot as plt

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries

# tstart = "2011-08-09 07:00"
# tend = "2011-08-09 10:00"
tstart = "2014-09-10 17:00"
tend = "2014-09-10 23:00"
tr = a.Time(tstart, tend)

results = Fido.search(tr, a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s") | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == 'SWPC'))

files = Fido.fetch(results)
goes = TimeSeries(files)

hek_results = results['hek']
flares_hek = hek_results[0]

fig, ax = plt.subplots()
goes.plot(axes=ax)
ax.axvline(parse_time(flares_hek['event_peaktime']).datetime)
ax.axvspan(
    parse_time(flares_hek['event_starttime']).datetime,
    parse_time(flares_hek['event_endtime']).datetime,
    alpha=0.2, label=flares_hek['fl_goescls']
)
ax.legend(loc=2)
ax.set_yscale('log')
ax.set_xlim(tr.start.to_datetime(), tr.end.to_datetime())

plt.show()
