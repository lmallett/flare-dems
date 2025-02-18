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
path = "D://emcubes_110809//emcubes_np//"
pathno094 = "D://emcubes_110809//emcubes_mchannel//m94_np//"
pathno131 = "D://emcubes_110809//emcubes_mchannel//m131_np//"
pathno171 = "D://emcubes_110809//emcubes_mchannel//m171_np//"
pathno193 = "D://emcubes_110809//emcubes_mchannel//m193_np//"
pathno211 = "D://emcubes_110809//emcubes_mchannel//m211_np//"
pathno335 = "D://emcubes_110809//emcubes_mchannel//m335_np//"
tstart = "2011-08-09 07:00"
tend = "2011-08-09 10:00"
flarestart = parse_time("2011-08-09 07:48:15.57").datetime
flareend = parse_time("2011-08-09 09:29:51.57").datetime
flaredate = "2011-08-09"
xregion = [410,601]
yregion = [410,611]

# path = "D://emcubes_140910//emcubes_np//"
# pathno335= "D://emcubes_140910//emcubes_no_335//"
# # xregion = [250,750]
# # yregion = [350,600]
# xregion = [375,575]
# yregion = [400,480]
# tstart = "2014-09-10 16:00"
# tend = "2014-09-10 23:00"
# flarestart = parse_time('2014-09-10 17:15:02.57').datetime
# flareend = parse_time('2014-09-10 18:44:38.57').datetime
# flaredate = "2014-09-10"

########################################################
########################################################

xstart = xregion[0]
xend = xregion[1] - 1
ystart = yregion[0]
yend = yregion[1] - 1

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

        elif file[-4:] == ".npz":
            vars = np.load(path + file)
            satmap = np.sum((vars["satmap"])[:, ystart:yend, xstart:xend], axis=0)
            statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]
            statuscube = np.ma.masked_where(
                satmap != 0, statuscube
            )  # exclude pixels that were saturated

        numfails[i] = np.count_nonzero(statuscube)

    return numfails

total_pxls = (yend - ystart) * (xend - xstart)
print(total_pxls)

def generate_fails():

    print("generating num_fails")
    num_fails = numfails(path)
    # num_fails_no094 = numfails(pathno094)
    # num_fails_no131 = numfails(pathno131)
    # num_fails_no171 = numfails(pathno171)
    # num_fails_no193 = numfails(pathno193)
    # num_fails_no211 = numfails(pathno211)
    num_fails_no335 = numfails(pathno335)
    np.savez("num_fails 2011", num_fails=num_fails, num_fails_no335=num_fails_no335)
    # np.savez("num_fails", num_fails = num_fails,
    #         num_fails_no094 = num_fails_no094,
    #         num_fails_no131 = num_fails_no131,
    #         num_fails_no171 = num_fails_no171,
    #         num_fails_no193 = num_fails_no193,
    #         num_fails_no211 = num_fails_no211,
    #         num_fails_no335 = num_fails_no335,
    #          )

def load_fails():
    global num_fails
    global num_fails_no335

    print("reading num_fails")
    vars = np.load("num_fails 2011.npz")
    num_fails = vars["num_fails"]
    # num_fails_no094 = vars["num_fails_no094"]
    # num_fails_no131 = vars["num_fails_no131"]
    # num_fails_no171 = vars["num_fails_no171"]
    # num_fails_no193 = vars["num_fails_no193"]
    # num_fails_no211 = vars["num_fails_no211"]
    num_fails_no335 = vars["num_fails_no335"]

########################################################
########################################################
# This part of the code creates GOES flux plots for a solar flare as desired.
# Mostly copied from the following documentation:
# https://docs.sunpy.org/en/latest/generated/gallery/time_series/goes_hek_m25.html

filelist = os.listdir(path)

def plot_two():

    fig, axs = plt.subplots(
        2, 1,
        sharex  = True, 
        figsize = (14, 8), 
        gridspec_kw = {"height_ratios": [2, 1]}
    )

    tr = a.Time(flarestart, flareend)
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

    axs[0].set_xlim(tr.start.to_datetime(), tr.end.to_datetime())
    axs[0].set_yticks(ticks=[1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2], minor=True)
    axs[0].set_yscale("log")
    axs[0].yaxis.grid(True)

    axs[0].set_xlabel("Time (HH:MM)")
    axs[0].set_ylabel(
        r"SXR Flux $[\mathrm{W} \; \mathrm{m}^{-2}]$", labelpad=80, # rotation="horizontal",
    )

    axs[0].legend(labels = ["0.5--4.0 \AA","1.0--8.0 \AA", "Excluded from analysis due to saturation"])

    goesclass = axs[0].twinx()
    goesclass.set_yticks([])
    # goesclass.set_ylabel("GOES \nFlare \nClass", rotation="horizontal", va="center", ha="left", labelpad=40)

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
        "With all channels": num_fails/total_pxls,
        # "Without 94": num_fails_no094,
        # "Without 131": num_fails_no131,
        # "Without 171": num_fails_no171,
        # "Without 193": num_fails_no193,
        # "Without 211": num_fails_no211,
        "Without 335 \AA": num_fails_no335/total_pxls,
    }

    times_df = DataFrame(data, index=filelist_times)
    frametimes = TimeSeries(times_df)

    axs[1].set_prop_cycle(
        color=["#000000", 
            'm'],
            #    "#a666d2", "#979ce1", "#78cef0", "#ffc500", "#ff8400", "#ff0000"],
        # linestyle=["-", ":", "-", "--", ":", "--", "-"],
    )

    frametimes.plot(axes=axs[1], linewidth=1.5)

    axs[1].set_ylabel(
        "Num. Failed Pixels /\n Total Pixels", rotation="horizontal", va="center", labelpad=80
    )
    axs[1].set(ylim = [0,.6])
    axs[1].grid()
    axs[1].legend()

    plt.tight_layout()
    plt.show()

def plot_dem_fails():
    
    fig, axs = plt.subplots(figsize = (14,4))

    filelist_short = [(word[11:-4]).replace("_", ":") for word in filelist]
    filelist_times = [parse_time(t).datetime for t in filelist_short]

    data = {
        "With all channels": num_fails/total_pxls,
        # "Without 94": num_fails_no094,
        # "Without 131": num_fails_no131,
        # "Without 171": num_fails_no171,
        # "Without 193": num_fails_no193,
        # "Without 211": num_fails_no211,
        "Without 335 \AA": num_fails_no335/total_pxls,
    }

    times_df = DataFrame(data, index = filelist_times)
    frametimes = TimeSeries(times_df)

    axs.set_prop_cycle(color=["#000000", 'm'],)
    frametimes.plot(axes=axs, 
                    linewidth=1.5,
                    )
    
    axs.set_xlabel("Time (HH:MM)")
    axs.set_ylabel(
        "Num. Failed Pixels /\n Total Pixels",
        va="center",
        labelpad=25,
        #rotation="horizontal"
    )

    tr = a.Time(flarestart, flareend)
    axs.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))  # This just removes the axis text annotation giving the date
    axs.set(xlim = [tr.start.to_datetime(), tr.end.to_datetime()], ylim = [0, .6])
    axs.grid()

    axs.legend(title = flaredate)
    plt.tight_layout()
    plt.show()

# generate_fails()
load_fails()
plot_dem_fails()