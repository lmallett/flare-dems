# (c) Lucien Mallett
# This program produces a GOES flux plot over time
# As well as showing how many DEM inversion failures fail over time

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

mpl.use("TkAgg")
plt.rcParams["text.usetex"] = True

plt.rc("text", usetex=True)
font = {"family": "serif", "size": 16}  # , 'serif': ['computer modern roman']}
plt.rc("font", **font)
# plt.rc('legend',**{'fontsize':14})
# plt.rc('text.latex', preamble = r'\usepackage[utf8]{inputenc}')

from cycler import cycler

##############################################
##############################################
# Path to folder containing .npz files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

# Path to folder containing .sav files, which contain image data at different timesteps
# path = "D:\\emcubes_110809\\emcubes\\"
# path = "D:\\emcubes_no_335_110809\\emcubes_no_335\\"
path = "D://emcubes_110809//emcubes_np//"
pathno094 = "D://emcubes_110809//emcubes_mchannel//m94_np//"
pathno131 = "D://emcubes_110809//emcubes_mchannel//m131_np//"
pathno171 = "D://emcubes_110809//emcubes_mchannel//m171_np//"
pathno193 = "D://emcubes_110809//emcubes_mchannel//m193_np//"
pathno211 = "D://emcubes_110809//emcubes_mchannel//m211_np//"
pathno335 = "D://emcubes_110809//emcubes_mchannel//m335_np//"

xregion = [350, 650]
yregion = [275, 775]

# path = "D:\\emcubes_140910\\emcubes\\"
# path = "D:\\emcubes_140910\\emcubes_no_335\\"
# xregion = [250,750]
# yregion = [350,750]

########################################################
########################################################
# Variables
########################################################
########################################################

xstart = xregion[0]
xend = xregion[1] - 1
ystart = yregion[0]
yend = yregion[1] - 1

tstart = "2011-08-09 07:00"
tend = "2011-08-09 10:00"
diffractionstart = parse_time("2011-08-09 08:01:39.57Z").datetime
diffractionend = parse_time("2011-08-09 08:06:18.01Z").datetime
flarestart = parse_time("2011-08-09 07:48:15.57").datetime
flareend = parse_time("2011-08-09 09:29:51.57").datetime

# tstart = "2014-09-10 16:00"
# tend = "2014-09-10 23:00"
# diffractionstart = parse_time('2014-09-10 17:30:26.57Z').datetime
# diffractionend = parse_time('2014-09-10 17:32:14.57Z').datetime
# flarestart = parse_time('2014-09-10 17:15:02.57').datetime
# flareend = parse_time('2014-09-10 18:44:38.57').datetime

##############################################
##############################################
# Data
##############################################
##############################################


def numfails(path):

    filelist = os.listdir(path)
    numfails = np.zeros_like(filelist, dtype=int)
    for i, file in enumerate(filelist):
        print(i)

        if file[-4:] == ".sav":
            vars = sp.io.readsav(path + file)
            satmap = np.sum((vars.satmap)[:, ystart:yend, xstart:xend], axis=0)
            statuscube = (vars.statuscube)[ystart:yend, xstart:xend]
            statuscube = np.ma.masked_where(
                satmap != 0, statuscube
            )  # exclude pixels that were saturated

            # emcube = np.sum(vars.emcube[:,ystart:yend, xstart:xend], axis = 0)              # creates 2D emcube array
            # emcube = emcube + np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis = 0)
            # numfails[i] = np.count_nonzero(emcube == 0) # counterintuitively, this counts zeros.
        elif file[-4:] == ".npz":
            vars = np.load(path + file)
            satmap = np.sum((vars["satmap"])[:, ystart:yend, xstart:xend], axis=0)
            statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]
            statuscube = np.ma.masked_where(
                satmap != 0, statuscube
            )  # exclude pixels that were saturated

        numfails[i] = np.count_nonzero(statuscube)

    return numfails


try:
    vars = np.load("num_fails.npz")
    num_fails = vars["num_fails"]
    num_fails_no094 = vars["num_fails_no094"]
    num_fails_no131 = vars["num_fails_no131"]
    num_fails_no171 = vars["num_fails_no171"]
    num_fails_no193 = vars["num_fails_no193"]
    num_fails_no211 = vars["num_fails_no211"]
    num_fails_no335 = vars["num_fails_no335"]

except:
    num_fails = numfails(path)
    num_fails_no094 = numfails(pathno094)
    num_fails_no131 = numfails(pathno131)
    num_fails_no171 = numfails(pathno171)
    num_fails_no193 = numfails(pathno193)
    num_fails_no211 = numfails(pathno211)
    num_fails_no335 = numfails(pathno335)
    np.savez("num_fails", num_fails=num_fails, num_fails_no335=num_fails_no335)

# num_fails = numfails(path)
# num_fails_no094 = numfails(pathno094)
# num_fails_no131 = numfails(pathno131)
# num_fails_no171 = numfails(pathno171)
# num_fails_no193 = numfails(pathno193)
# num_fails_no211 = numfails(pathno211)
# num_fails_no335 = numfails(pathno335)
# np.savez("num_fails", num_fails = num_fails,
#         num_fails_no094 = num_fails_no094,
#         num_fails_no131 = num_fails_no131,
#         num_fails_no171 = num_fails_no171,
#         num_fails_no193 = num_fails_no193,
#         num_fails_no211 = num_fails_no211,
#         num_fails_no335 = num_fails_no335,
#          )


filelist = os.listdir(path)

########################################################
########################################################
# Plotting
########################################################
########################################################

fig, axs = plt.subplots(
    2, 1, sharex=True, figsize=(14, 8), gridspec_kw={"height_ratios": [2, 1]}
)
# fig.subplots_adjust(hspace = 0)

########################################################
# This part of the code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

tr = a.Time(tstart, tend)
results = Fido.search(
    tr,
    a.Instrument("XRS") & a.goes.SatelliteNumber(15) & a.Resolution("flx1s")
    | a.hek.FL & (a.hek.FL.GOESCls > "X1.0") & (a.hek.FRM.Name == "SWPC"),
)
files = Fido.fetch(results)
goes = TimeSeries(files)
plot = goes.plot(columns=["xrsa", "xrsb"], axes=axs[0])

hek_results = results["hek"]
flares_hek = hek_results[0]
axs[0].axvline(parse_time(flares_hek["event_peaktime"]).datetime)
axs[0].axvspan(
    flarestart,
    flareend,
    alpha=0.2,
    label=np.str_("Data Gathered"),
    # label= flares_hek['fl_goescls']
)

# Show red region where we are excluding data
axs[0].axvspan(
    diffractionstart,
    diffractionend,
    alpha=0.2,
    facecolor="r",
    label=np.str_("Data Excluded"),
)

axs[0].legend(loc=1)

axs[0].set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
axs[0].set_yticks(ticks=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2], minor=True)
axs[0].set_yscale("log")
axs[0].yaxis.grid(True)

axs[0].set_xlabel("Time")
axs[0].set_ylabel(
    r"SXR Flux $[\mathrm{W} \; \mathrm{m}^{-2}]$", rotation="horizontal", labelpad=80
)

goesclass = axs[0].twinx()
goesclass.set_ylabel(
    "GOES \nFlare \nClass", rotation="horizontal", va="center", ha="left", labelpad=40
)
goesclass.set_yticks([])

# # the GOES class labels are contained in the class XRSTimeSeries(GenericTimeSeries) in sunpy so we need to rewrite this
# # section to make the flare classes are on the left side if you want
# for t in plot.texts:
#     t.set_visible(False)
# labels = ['A', 'B', 'C', 'M', 'X']
# centers = np.logspace(-7.5, -3.5, len(labels))
# for value, label in zip(centers, labels):
#     axs[0].text(-0.02, value, label, transform=axs[0].get_yaxis_transform(), horizontalalignment='center')

########################################################
# Fails of DEM

filelist_short = [(word[11:-4]).replace("_", ":") for word in filelist]
filelist_times = [parse_time(t).datetime for t in filelist_short]

# data = {"With 335": fails2011, "Without 335": fails2011_no335, "With 335 No Sat": fails2011_nosat, "Without 335 No Sat": fails2011_no335_nosat}
data = {
    "With all channels": num_fails,
    "Without 94": num_fails_no094,
    "Without 131": num_fails_no131,
    "Without 171": num_fails_no171,
    "Without 193": num_fails_no193,
    "Without 211": num_fails_no211,
    "Without 335": num_fails_no335,
}

times_df = DataFrame(data, index=filelist_times)
frametimes = TimeSeries(times_df)


axs[1].set_prop_cycle(
    color=["#000000", "#a666d2", "#979ce1", "#78cef0", "#ffc500", "#ff8400", "#ff0000"],
    linestyle=["-", ":", "-", "--", ":", "--", "-"],
)

frametimes.plot(axes=axs[1], linewidth=1.5)
axs[1].legend(ncol = 2)
axs[1].set_ylabel(
    "Num. DEM \nInversion Failures", rotation="horizontal", va="center", labelpad=80
)
axs[1].grid()


plt.tight_layout()
plt.show()

# print(type(flares_hek['fl_goescls']))
# print(flares_hek['fl_goescls'])
# print(type(flares_hek['fl_goescls']))
# print(np.str_("Data Used"))
# print(type(np.str_("Data Used")))
