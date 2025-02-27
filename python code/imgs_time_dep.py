# (c) Lucien Mallett 2024

import scipy as sp
import numpy as np
import os

import matplotlib
# matplotlib.use("Agg")
matplotlib.use('TkAgg')

# from matplotlib import rc
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy import visualization

##############################################
##############################################
# Path to folder containing .sav files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

# pathresp = "aia_resp_110809.sav"
# path = "D://emcubes_110809//emcubes_np//"
# # path = "D://emcubes_no_335_110809//emcubes_no_335_np//"
# # dir_emcubes = savedir + "emcubes_no335//"
# # xregion = [350,651]
# # yregion = [305,746]
# # xregion = [350,650]
# # yregion = [275,775]
# xregion = [410,601]
# yregion = [410,611]

# savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2011//"

##############################################

pathresp = "aia_resp_140910.sav"
path = "D://emcubes_140910//emcubes_np//"
xregion = [250,750]
yregion = [350,600]
savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2014//"

##############################################
##############################################
# Initialization
##############################################
##############################################

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

# load in response functions as 20x1 arrays
respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

logte= (respinfo["aia_tresp"][0])[0]
r94  = (respinfo["aia_tresp"][0])[2]
r131 = (respinfo["aia_tresp"][0])[3]
r171 = (respinfo["aia_tresp"][0])[4]
r193 = (respinfo["aia_tresp"][0])[5]
r211 = (respinfo["aia_tresp"][0])[6]
r335 = (respinfo["aia_tresp"][0])[7]

allwaves  = [94,131,171,193,211,335]
filelist = os.listdir(path)

##############################################
##############################################
# Plot
##############################################
##############################################

def exptimes(wavelength = allwaves):
    """Produces a graph of the exposure times for each channel over time.

    Keyword Arguments:
        wavelength -- An integer or list of integers of channel wavelengths in Angstroms. (default: {allwaves}, all available channels.)
    """

    if type(wavelength) != list:
        wavelength = [wavelength]

    n = len(filelist)

    exptimes = np.zeros([n, 6])
    times = []
    print("done")

    for i in range(n):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        time = file[11:-4].replace("_",":")
        times.append(time)

        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)
        exptimes[i] = np.copy(vars.exptimeout)
        print(filename)

    ts = TimeSeries(time = times)
    ts.sort('time')
    ts['times_94'] = exptimes[:,0]
    ts['times_131'] = exptimes[:,1]
    ts['times_171'] = exptimes[:,2]
    ts['times_193'] = exptimes[:,3]
    ts['times_211'] = exptimes[:,4]
    ts['times_335'] = exptimes[:,5]

    print(ts)

    print(ts.time.format)
    fig, ax = plt.subplots()

    with visualization.time_support(format = 'isot', simplify=True):  

        ax.plot(ts.time, ts['times_94'], label = '94')
        ax.plot(ts.time, ts['times_131'], label = '131')
        ax.plot(ts.time, ts['times_171'], label = '171')
        ax.plot(ts.time, ts['times_193'], label = '193')
        ax.plot(ts.time, ts['times_211'], label = '211')
        ax.plot(ts.time, ts['times_335'], label = '335')

        # ax.xlabel('Julian Date')
        # ax.ylabel('SAP Flux (e-/s)')
        ax.grid()
    plt.legend()
    
    plt.show()

######################################################################
######################################################################

def num_sat():
    """Produces a "number saturated vs. time" graph for a given data set.
    """

    n = len(filelist)

    numSat = np.zeros([n, 6])
    norms = np.zeros([n, 6])

    times = []

    for i in range(n):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        time = file[11:-4].replace("_",":")
        times.append(time)

        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)
        sums = np.sum(vars.satmap, axis = (1,2))

        numSat[i] = sums

        # Catch for when there are no saturated channels by setting them to 0
        norms[i] = np.nan_to_num(sums / np.linalg.norm(sums))

        print(sums)

    ts = TimeSeries(time = times)
    ts['times_94'] = numSat[:,0]
    ts['times_131'] = numSat[:,1]
    ts['times_171'] = numSat[:,2]
    ts['times_193'] = numSat[:,3]
    ts['times_211'] = numSat[:,4]
    ts['times_335'] = numSat[:,5]

    ts['norm_94'] = norms[:,0]
    ts['norm_131'] = norms[:,1]
    ts['norm_171'] = norms[:,2]
    ts['norm_193'] = norms[:,3]
    ts['norm_211'] = norms[:,4]
    ts['norm_335'] = norms[:,5]

    ts.sort('time')
    print(ts)
    
    # print("Average times when there are saturated pixels:")
    # print("94A:" np.average())

    fig, ax = plt.subplots(2, sharex = True)

    with visualization.time_support(format = 'isot', simplify=True):  

        ax[0].plot(ts.time, ts['times_94'], label = "$94 \\mathrm{\AA}$", linestyle = 'dashed', c = 'r')
        ax[0].plot(ts.time, ts['times_131'], label = "$131 \\mathrm{\AA}$", c = 'y')
        ax[0].plot(ts.time, ts['times_171'], label = "$171 \\mathrm{\AA}$", linestyle = 'dashed', c = 'g')
        ax[0].plot(ts.time, ts['times_193'], label = "$193 \\mathrm{\AA}$", c = 'c')
        ax[0].plot(ts.time, ts['times_211'], label = "$211 \\mathrm{\AA}$", linestyle = 'dashed', c = 'm')
        ax[0].plot(ts.time, ts['times_335'], label = "$335 \\mathrm{\AA}$", c = 'k')
        ax[0].grid()
        ax[0].legend()

        ax[1].plot(ts.time, ts['norm_94'], label = "$94 \\mathrm{\AA}$", linestyle = 'dashed', c = 'r')
        ax[1].plot(ts.time, ts['norm_131'], label = "$131 \\mathrm{\AA}$", c = 'y')
        ax[1].plot(ts.time, ts['norm_171'], label = "$171 \\mathrm{\AA}$", linestyle = 'dashed', c = 'g')
        ax[1].plot(ts.time, ts['norm_193'], label = "$193 \\mathrm{\AA}$", c = 'c')
        ax[1].plot(ts.time, ts['norm_211'], label = "$211 \\mathrm{\AA}$", linestyle = 'dashed', c = 'm')
        ax[1].plot(ts.time, ts['norm_335'], label = "$335 \\mathrm{\AA}$", c = 'k')
        ax[1].grid()
        ax[1].legend()
        ax[1].set(ylim = [-0.1,1.1])

    plt.show()

def percent_dem_fails_saturated():

    # for each channel make a list
    n = len(filelist)
    times = []
    percent_fails = np.empty((6, n))
    print(percent_fails.shape)

    # for every frame, count the intersections between statuscube (2 and 10) and satmap
    for i in range(n):  # for each timestamp

        # read in data & modify it
        file = filelist[i]
        filename = path + file
        
        time = file[11:-4].replace("_",":")
        times.append(time)

        vars = np.load(path + file)
        satmap = ((vars["satmap"])[:, ystart:yend, xstart:xend]).astype(bool)   # makes T/F cube, length x width x 6 channels, to describe if a pixel is saturated in a channel

        statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]
        statuscube = np.isin(statuscube, [2,10])                                # makes T/F array, length x width, to describe if an error code is 2 or 10

        num_fail = np.sum(statuscube)   # number of failed pixels, scalar

        num_sat_fail = np.sum(np.sum(np.logical_and(satmap, statuscube), axis = 1), axis = 1) # number of pixels that failed AND are saturated, scalar

        print(num_fail)
        print(num_sat_fail)

        # then, this is a 6 x n array where each of the 6 is an AIA channel and n is over time
        percent_fails[:,i] = np.divide(num_sat_fail, num_fail, out=np.zeros_like(num_sat_fail, dtype=float), where = num_fail != 0)

        # print(percent_fails_row.shape)
        # print(sat_and_fail)
        # print(num_sat)
        # print(percent_fails[:,i])


    print(np.shape(times))
    print(np.shape(percent_fails[0,:]))

    ts = TimeSeries(time = times)
    ts['times_94'] = percent_fails[0,:]
    ts['times_131'] = percent_fails[1,:]
    ts['times_171'] = percent_fails[2,:]
    ts['times_193'] = percent_fails[3,:]
    ts['times_211'] = percent_fails[4,:]
    ts['times_335'] = percent_fails[5,:]

    ts.sort('time')
    ts.write(savedir + 'table.csv', format='ascii')

    print(ts)

    fig, ax = plt.subplots()

    with visualization.time_support(format = 'isot', simplify=True):  

        ax.plot(ts.time, ts['times_94'], label = "$94 \\mathrm{\AA}$", linestyle = 'dashed', c = 'r')
        ax.plot(ts.time, ts['times_131'], label = "$131 \\mathrm{\AA}$", c = 'y')
        ax.plot(ts.time, ts['times_171'], label = "$171 \\mathrm{\AA}$", linestyle = 'dashed', c = 'g')
        ax.plot(ts.time, ts['times_193'], label = "$193 \\mathrm{\AA}$", c = 'c')
        ax.plot(ts.time, ts['times_211'], label = "$211 \\mathrm{\AA}$", linestyle = 'dashed', c = 'm')
        ax.plot(ts.time, ts['times_335'], label = "$335 \\mathrm{\AA}$", c = 'k')

        ax.set(xlim =[ts.time[0], ts.time[-1]],
               ylim = [0,1],
               xlabel = "Time (HH:MM)",
               ylabel = "Proportion of failures that are saturated in a given channel")
        ax.grid()
        ax.legend()
        # ax.set_xlim(ts.time[0], ts.time[-1])
        # ax.set_xlabel("Time (HH:MM)")
        ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))  # This just removes the axis text annotation giving the date
    
    plt.show()

