# Lucien Mallett
# This code produces GOES flux plots. Requires user to manually save the result.
# Mostly copied from https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['text.usetex'] = True
matplotlib.use("TkAgg")

plt.rc('text',usetex = True)
font = {'family':'serif','size':16}
plt.rc('font',**font)

##############################################
##############################################
# (Un)comment or adjust these lines to switch which flare to plot the data for.

# # 2011 flare
# tstart = parse_time("2011-08-09 07:45").datetime
# tend = parse_time("2011-08-09 09:30").datetime
# impulsetime = parse_time("2011-08-09 08:01:39.57Z").datetime
# flaredate = "2011-08-09"

# 2014 flare
tstart = parse_time("2014-09-10 17:00").datetime
tend = parse_time("2014-09-10 19:00").datetime
impulsetime = parse_time("2014-09-10 17:22:00.00Z").datetime
flaredate = "2014-09-10"

##############################################
##############################################

tr = a.Time(tstart, tend)

results = Fido.search(tr, a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s") | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == 'SWPC'))

files = Fido.fetch(results)
goes = TimeSeries(files)

hek_results = results['hek']
flares_hek = hek_results[0]

fig, ax = plt.subplots(figsize = (12,5.25))
plot = goes.plot(axes = ax)

ax.set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
ax.set_xlabel("Time (HH:MM)")
ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))  # Removes the axis text annotation giving the date so we can put it in the legend

ax.set_yscale('log')
ax.set_ylabel(r'SXR Flux $[\mathrm{W} \; \mathrm{m}^{-2}]$')
ax.legend(title = flaredate, labels = ["0.5--4.0 \AA","1.0--8.0 \AA"], framealpha = 1)
ax.set_ylim(ymin= 1e-6, ymax = 1e-2)    # restricting y-axis to only show C-X class 

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

gradualtime = parse_time(flares_hek['event_peaktime']).datetime
ax.axvline(x = impulsetime)
ax.axvline(x = gradualtime)

ax.annotate(text = "Pre-flare \n phase", 
            xy = (tstart + (impulsetime - tstart)/2, 0.75), 
            xycoords = ('data', 'figure fraction'),
            ha = "center",
            # bbox = dict(boxstyle = 'square, pad = 1',
                    #   fc = 'red',
                    #   ec = 'none')
            )

##############################################

# (Un)comment or adjust lines to switch flares.

# 2011 flare
# ax.annotate(text = "Impulsive \n phase", 
#             xy = (impulsetime + (gradualtime - impulsetime)/3, 0.75), 
#             xycoords = ('data', 'figure fraction'),
#             xytext = (70, 0), 
#             textcoords = "offset points",
#             ha = 'center',
#             arrowprops=dict(arrowstyle="->", connectionstyle="arc3")
#             )

# # 2014 flare
ax.annotate(text = "Impulsive \n phase", 
            xy = (impulsetime + (gradualtime - impulsetime)/2, 0.75), 
            xycoords = ('data', 'figure fraction'),
            ha = 'center',
            )

##############################################

ax.annotate(text = "Gradual \n phase", 
            xy = (gradualtime + (parse_time(tend).datetime - gradualtime)/3, 0.75), 
            xycoords = ('data', 'figure fraction'),
            xytext = (10, 0), 
            textcoords = "offset points",
            ha = 'center'
            )

plt.show()