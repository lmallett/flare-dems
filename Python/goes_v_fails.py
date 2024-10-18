
import os
import matplotlib.pyplot as plt

import numpy as np
import scipy as sp


import matplotlib as mpl
mpl.use('TkAgg')
plt.rcParams['text.usetex'] = True

# This part of the code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

# plt.rcParams['text.usetex'] = True
plt.rc('text',usetex=True)
font = {'family':'serif','size':16}#, 'serif': ['computer modern roman']}
plt.rc('font',**font)
# plt.rc('legend',**{'fontsize':14})
# plt.rc('text.latex', preamble = r'\usepackage[utf8]{inputenc}')

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries


##############################################
##############################################
# CHANGE THESE LINES!
##############################################
##############################################

# Path to folder containing .sav files, which contain image data at different timesteps
# path = "E:\\2011_08_09\\emcubes\\all\\"
# path = "E:\\2014_09_10\\emcubes\\all\\"
# path = "E:\\2011_08_09.tar\\2011_08_09\\emcubes\\all\\"


# path = "E:\\emcubes_110809\\"
path = "E:\\emcubes_110809\\emcubes\\"


xregion = [350,650]
yregion = [275,775]

# path = "E:\\emcubes_140910\\"
# xregion = [250,750]
# yregion = [350,750]   

##############################################
##############################################

filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1


##############################################
##############################################

# numfails = np.zeros_like(filelist)

# for i, file in enumerate(filelist):
#     vars = sp.io.readsav(path + file)
#     emcube = np.sum(vars.emcube[:,ystart:yend, xstart:xend], axis = 0) # creates 2D emcube array
#     numfails[i] = np.count_nonzero(emcube == 0) # counterintuitively, this counts zeros.

# print(numfails)

fails2011 = [649, 651, 693, 686, 691, 723, 614, 604, 592, 605, 536, 557,
    570, 571, 620, 610, 599, 613, 558, 594, 669, 587, 227, 115,
    277, 8094, 29386, 50155, 58047, 66908, 75908, 75239, 49020,
    58687, 54451, 48449, 29030, 28093, 34331, 26049, 17225, 13009,
    20963, 19380, 7789, 18449, 13152, 6869, 19792, 11071, 4042,
    2834, 7968, 5465, 1665, 8194, 2122, 1596, 4623, 2401, 1145,
    4146, 987, 1022, 2372, 947, 956, 856, 1000, 1107, 980, 860,
    1057, 1045, 648, 615, 780, 855, 982, 537, 822, 1023, 823,
    1040, 1185, 1627, 1163, 1106, 1313, 1138, 1206, 1051, 1058,
    1218, 1049, 912, 1204, 969, 968, 895, 922, 1122, 886, 782,
    1184, 854, 676, 617, 1365, 1409, 1019, 1159, 1266, 961, 1347,
    1070, 991, 1251, 1105]


########################################################
# This part of the code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

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

########################################################

filelist_short = [word[11:-4] for word in filelist]
print(filelist_short)

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

fig, ax = plt.subplots(figsize = (12,5.25))

goes.plot(axes = ax)

ax.axvline(parse_time(flares_hek['event_peaktime']).datetime)
ax.axvspan(
    flarestart,
    flareend,
    alpha=0.2,
    label = np.str_("Data Gathered")
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

# ax.legend()
ax.legend(loc=1)
ax.set_yscale('log')
ax.set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
ax.set_ylabel(r'SXR Flux $(\mathrm{W} \; \mathrm{m}^{-2})$')
ax.set_xlabel("Time")
plt.tight_layout()
plt.show()
