import os
import numpy as np
import scipy as sp

from sunpy.net import Fido
from sunpy.net import attrs as a
from sunpy.time import parse_time
from sunpy.timeseries import TimeSeries
from pandas import DataFrame

import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('TkAgg')
plt.rcParams['text.usetex'] = True

# plt.rcParams['text.usetex'] = True
plt.rc('text',usetex=True)
font = {'family':'serif','size':16}#, 'serif': ['computer modern roman']}
plt.rc('font',**font)
# plt.rc('legend',**{'fontsize':14})
# plt.rc('text.latex', preamble = r'\usepackage[utf8]{inputenc}')

##############################################
##############################################
# Path strings and region selection
##############################################
##############################################

# Path to folder containing .sav files, which contain image data at different timesteps
# path = "E:\\emcubes_110809\\emcubes\\"
path = "E:\\emcubes_no_335_110809\\emcubes_no_335\\"

xregion = [350,650]
yregion = [275,775]

# path = "E:\\emcubes_140910\\emcubes\\"
# path = "E:\\emcubes_140910\\emcubes_no_335\\"
# xregion = [250,750]
# yregion = [350,750]   

##############################################
##############################################
# Data
##############################################
##############################################

# # Uncomment these lines to get results; save them as a list

# numfails = np.zeros_like(filelist)
# for i, file in enumerate(filelist):
#     print(i)
#     vars = sp.io.readsav(path + file)
#     emcube = np.sum(vars.emcube[:,ystart:yend, xstart:xend], axis = 0)              # creates 2D emcube array
#     # add in saturation map so that we don't include pixels that are saturated. since satmap is 1/0, excludes pixels that were saturated in DEM fail from count
#     emcube = emcube + np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis = 0)     

#     numfails[i] = np.count_nonzero(emcube == 0) # counterintuitively, this counts zeros.

# print(numfails)

##############################################
# This data is produced from numfails

fails2011 = [
    649, 651, 693, 686, 691, 
    723, 614, 604, 592, 605, 
    536, 557, 570, 571, 620, 
    610, 599, 613, 558, 594, 
    669, 587, 227, 115, 277, 
    8094, 29386, 50155, 58047, 66908, 
    75908, 75239, 49020, 58687, 54451, 
    48449, 29030, 28093, 34331, 26049, 
    17225, 13009, 20963, 19380, 7789, 
    18449, 13152, 6869, 19792, 11071, 
    4042, 2834, 7968, 5465, 1665, 
    8194, 2122, 1596, 4623, 2401, 
    1145, 4146, 987, 1022, 2372, 
    947, 956, 856, 1000, 1107, 
    980, 860, 1057, 1045, 648, 
    615, 780, 855, 982, 537, 
    822, 1023, 823, 1040, 1185, 
    1627, 1163, 1106, 1313, 1138, 
    1206, 1051, 1058, 1218, 1049, 
    912, 1204, 969, 968, 895, 
    922, 1122, 886, 782, 1184, 
    854, 676, 617, 1365, 1409, 
    1019, 1159, 1266, 961, 1347, 
    1070, 991, 1251, 1105]

fails2011_nosat = [
    649, 651, 693, 686, 691, 
    723, 614, 604, 592, 605, 
    536, 557, 570, 571, 620,
    610, 599, 613, 558, 594, 
    669, 587, 226, 110, 226, 
    6137, 28980, 49695, 57488, 66066, 
    75534, 75028, 48625,  58483, 54132, 
    48305, 28738, 27832, 34190, 25980, 
    17078, 12808,  20852, 19297, 7642, 
    18308, 13033, 6765, 19619, 11015, 
    3969,  2759, 7925, 5401, 1579, 
    8163, 2068, 1523, 4597, 2370, 
    1088, 4137, 952, 992, 2363, 
    905, 895, 802, 921, 977, 
    812, 822, 962, 939, 616, 
    606, 752, 813, 943, 530, 
    803, 1017, 819, 1039, 1183, 
    1625, 1163, 1105, 1311, 1137, 
    1206, 1051, 1056, 1217, 1049, 
    911, 1198, 969, 965, 895, 
    921, 1116, 886, 782, 1184, 
    854, 676, 615, 1364, 1408, 
    1018, 1159, 1266, 961, 1347, 
    1070, 991, 1251, 1105]

fails2011_no335 = [
    0, 0, 0, 2, 0, 
    1, 0, 0, 0, 0,
    0, 0, 1, 0, 0, 
    0, 1, 2, 0, 0,
    137, 10, 75, 72, 205, 
    7679, 27766, 44914, 47134, 50512,
    48476, 52186, 34035, 47735, 37446, 
    34646, 23780, 25338, 26136, 23037,
    16029, 11170, 20148, 17923, 7434, 
    18202, 12264, 6487, 18001, 10826,
    3784, 2246, 7990, 4827, 1186, 
    7996, 1933, 1032, 4676, 2285,
    614, 4197, 853, 565, 2245, 
    704, 534, 368, 515, 154,
    488, 478, 189, 521, 303, 
    248, 392, 196, 414, 413,
    152, 522, 517, 418, 655, 
    190, 453, 659, 4, 546,
    763, 381, 560, 13, 418, 
    493, 11, 454, 539, 348,
    503, 14, 318, 360, 4, 
    366, 288, 79, 171, 63,
    2, 50, 46, 0, 177, 
    9, 1, 81, 10]

fails2011_no335_nosat = [
    0, 0, 0, 2, 0,
    1, 0, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 1, 2, 0, 0,
    137, 10, 75, 72, 169,
    5841, 27374, 44443, 46590, 49661,
    48109, 51974, 33638, 47503, 37144,
    34483, 23470, 25113, 25992, 22959,
    15929, 11037, 20094, 17889, 7363,
    18117, 12174, 6429, 17940, 10802, 
    3743, 2201, 7963, 4788, 1147, 
    7989, 1918, 1010, 4667, 2273, 
    601, 4190, 840, 553, 2234, 
    697, 520, 353, 491, 130, 
    457, 477, 162, 499, 297, 
    239, 385, 161, 385, 405, 
    132, 518, 503, 418, 654, 
    185, 452, 657, 4, 546, 
    763, 380, 560, 7, 418, 
    493, 5, 454, 539, 348, 
    503, 10, 318, 360, 4, 
    366, 288, 78, 171, 63, 
    2, 50, 46, 0, 177, 
    9, 1, 81, 10]

########################################################
########################################################
# Variables
########################################################
########################################################

filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

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
########################################################
# Plotting
########################################################
########################################################

fig, axs = plt.subplots(2, 1, 
                        sharex = True, 
                        figsize = (15,7),
                        gridspec_kw= {
                            "height_ratios": [2,1]
                        }
                        )
# fig.subplots_adjust(hspace = 0)

########################################################
# This part of the code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

tr = a.Time(tstart, tend)
results = Fido.search(tr, a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s") | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == 'SWPC'))
files = Fido.fetch(results)
goes = TimeSeries(files)
plot = goes.plot(columns = ["xrsa", "xrsb"], axes = axs[0])

hek_results = results['hek']
flares_hek = hek_results[0]
axs[0].axvline(parse_time(flares_hek['event_peaktime']).datetime)
axs[0].axvspan(
    flarestart,
    flareend,
    alpha=0.2,
    label = np.str_("Data Gathered")
    # label= flares_hek['fl_goescls']
)
# Show red region where we are excluding data
axs[0].axvspan(
    diffractionstart,
    diffractionend,
    alpha=0.2,
    facecolor = 'r', 
    label = np.str_("Data Excluded")
)

axs[0].set_yticks(ticks = [1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2], minor = True)
axs[0].yaxis.grid(True)

axs[0].set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
axs[0].set_ylabel(r'SXR Flux $[\mathrm{W} \; \mathrm{m}^{-2}]$', rotation = 'horizontal', labelpad = 80)
axs[0].legend(loc = 1)

axs[0].set_yscale('log')
axs[0].set_xlabel("Time")

goesclass = axs[0].twinx()
goesclass.set_ylabel("GOES \nFlare \nClass", rotation = 'horizontal', va = 'center', ha = 'left', labelpad = 40)
goesclass.set_yticks([])

# # the GOES class labels are contained in the class XRSTimeSeries(GenericTimeSeries) in sunpy so we need to rewrite this
# # section to make the flare classes are on the left side if you want
# for t in plot.texts:
#     t.set_visible(False)
# labels = ['A', 'B', 'C', 'M', 'X']
# centers = np.logspace(-7.5, -3.5, len(labels))
# for value, label in zip(centers, labels):
#     axs[0].text(-0.02, value, label, transform=axs[0].get_yaxis_transform(), horizontalalignment='center')

for t in plot.texts:
    t.set_visible(False)
labels = ['A', 'B', 'C', 'M', 'X']
centers = np.logspace(-7.5, -3.5, len(labels))
for value, label in zip(centers, labels):
    axs[0].text(1.02, value, label, transform=axs[0].get_yaxis_transform(), horizontalalignment='center', verticalalignment = 'top')

########################################################
# Fails of DEM

filelist_short = [(word[11:-4]).replace('_', ':') for word in filelist]
filelist_times = [parse_time(t).datetime for t in filelist_short]

data = {"With 335": fails2011, "Without 335": fails2011_no335, "With 335 No Sat": fails2011_nosat, "Without 335 No Sat": fails2011_no335_nosat}
data = {"With 335": fails2011, "Without 335": fails2011_no335}

times_df = DataFrame(data, index = filelist_times)
frametimes = TimeSeries(times_df)

frametimes.plot(axes = axs[1])
axs[1].legend()
axs[1].set_ylabel("Num. DEM \nInversion Failures", rotation = 'horizontal', va = 'center', labelpad = 80)

plt.tight_layout()
plt.show()

# print(type(flares_hek['fl_goescls']))
# print(flares_hek['fl_goescls'])
# print(type(flares_hek['fl_goescls']))
# print(np.str_("Data Used"))
# print(type(np.str_("Data Used")))