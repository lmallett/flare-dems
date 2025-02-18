# This code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
matplotlib.use("TkAgg")


plt.rc('text',usetex=True)
font = {'family':'serif','size':16}#, 'serif': ['computer modern roman']}
plt.rc('font',**font)
# plt.rc('legend',**{'fontsize':14})
# plt.rc('text.latex', preamble = r'\usepackage[utf8]{inputenc}')

##############################################
##############################################

tstart = "2011-08-09 07:00"
tend = "2011-08-09 10:00"
diffractionstart = parse_time('2011-08-09 08:01:39.57Z').datetime
diffractionend = parse_time('2011-08-09 08:06:18.01Z').datetime
flarestart = parse_time('2011-08-09 07:48:15.57').datetime
flareend = parse_time('2011-08-09 09:29:51.57').datetime
flaredate = "2011-08-09"

# tstart = "2014-09-10 16:00"
# tend = "2014-09-10 23:00"
# diffractionstart = parse_time('2014-09-10 17:30:26.57Z').datetime
# diffractionend = parse_time('2014-09-10 17:32:14.57Z').datetime
# flarestart = parse_time('2014-09-10 17:15:02.57').datetime
# flareend = parse_time('2014-09-10 18:44:38.57').datetime
# flaredate = "2014-09-10"

##############################################
##############################################

# tr = a.Time(tstart, tend)
tr = a.Time(flarestart, flareend)

results = Fido.search(tr, a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s") | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == 'SWPC'))

files = Fido.fetch(results)
goes = TimeSeries(files)

hek_results = results['hek']
flares_hek = hek_results[0]

fig, ax = plt.subplots(figsize = (12,5.25))
plot = goes.plot(axes = ax)

ax.set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
ax.set_xlabel("Time (HH:MM)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))  # This just removes the axis text annotation giving the date

ax.set_yscale('log')
ax.set_ylabel(r'SXR Flux $[\mathrm{W} \; \mathrm{m}^{-2}]$')
ax.legend(title = flaredate, labels = ["0.5--4.0 \AA","1.0--8.0 \AA", "Excluded from analysis due to saturation"])
ax.set_ylim(ymin= 1e-6, ymax = 1e-2)

goesclass = ax.twinx()
goesclass.set_yticks([])

# the GOES class labels are contained in the class XRSTimeSeries(GenericTimeSeries) in sunpy so we need to rewrite this
# section to make the flare classes adjust with range
for t in plot.texts:
    t.set_visible(False)
labels = ['C', 'M', 'X']
centers = np.logspace(-5.5, -3.5, len(labels))
for value, label in zip(centers, labels):
    goesclass.text(1.02, value, label, transform= plot.get_yaxis_transform(), horizontalalignment='center')

plt.tight_layout()
plt.show()