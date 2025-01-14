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

# Path to response functions
pathresp = "aia_resp_full.sav"

path = "D:\\emcubes_110809\\emcubes\\"
# path = "D:\\emcubes_no_335_110809\\emcubes_no_335\\"
xregion = [350,650]
yregion = [275,775]

# path = "D:\\emcubes_140910\\emcubes\\"
# path = "D:\\emcubes_140910\\emcubes_no_335\\"
# xregion = [250,750]
# yregion = [350,750]

savedir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2014//"

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

logte = (((respinfo['r'])['A94'])[0])['logte'][0]
r94 = (((respinfo['r'])['A94'])[0])['tresp'][0]
r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]

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
