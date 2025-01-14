# This code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries

import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True

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

# tstart = "2014-09-10 16:00"
# tend = "2014-09-10 23:00"
# diffractionstart = parse_time('2014-09-10 17:30:26.57Z').datetime
# diffractionend = parse_time('2014-09-10 17:32:14.57Z').datetime
# flarestart = parse_time('2014-09-10 17:15:02.57').datetime
# flareend = parse_time('2014-09-10 18:44:38.57').datetime

##############################################
##############################################

tr = a.Time(tstart, tend)

results = Fido.search(tr, a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s") | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == 'SWPC'))

files = Fido.fetch(results)
goes = TimeSeries(files)

hek_results = results['hek']
flares_hek = hek_results[0]

fig, ax = plt.subplots(figsize = (12,5.25))
goes.plot(axes = ax)

ax.legend(loc = 1)

ax.set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
ax.set_yscale('log')

ax.set_xlabel("Time")
ax.set_ylabel(r'SXR Flux $[\mathrm{W} \; \mathrm{m}^{-2}]$', rotation = 'horizontal', labelpad = 80)

ax.axvline(parse_time(flares_hek['event_peaktime']).datetime)
ax.axvspan(
    flarestart,
    flareend,
    alpha=0.2,
    label = np.str_("Data Gathered")
    # label= flares_hek['fl_goescls']
)

# Show red region where we are excluding data
ax.axvspan(
    diffractionstart,
    diffractionend,
    alpha=0.2,
    facecolor = 'r',
    label = np.str_("Data Excluded")
)

goesclass = ax.twinx()
goesclass.set_ylabel("GOES \nFlare \nClass", rotation = 'horizontal', va = 'center', ha = 'left', labelpad = 40)
goesclass.set_yticks([])

plt.tight_layout()
plt.show()

# print(type(flares_hek['fl_goescls']))
# print(flares_hek['fl_goescls'])
# print(type(flares_hek['fl_goescls']))
# print(np.str_("Data Used"))
# print(type(np.str_("Data Used")))