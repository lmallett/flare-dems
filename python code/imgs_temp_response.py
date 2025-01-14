#  (c) Lucien Mallett 2024

import scipy as sp
import numpy as np
import os

import matplotlib
# matplotlib.use("Agg")
matplotlib.use('TkAgg')

matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import matplotlib.pyplot as plt

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
# Plotting
##############################################
##############################################

def norm_temp_responses():
    """Plot the temperature response functions for AIA, normalized, to show where each channel is most sensitive."""

    # load in response functions as 101x1 arrays
    respinfo = sp.io.readsav("aia_response_2011.sav", python_dict = False, verbose = True)
    r94 = respinfo['r94']
    r131 = respinfo['r131']
    r171 = respinfo['r171']
    r193 = respinfo['r193']
    r211 = respinfo['r211']
    r335 = respinfo['r335']
    logte = ((respinfo['r'])['logte'])[0]

    #normalize the functions
    r94 = r94 / max(r94)
    r131 = r131 / max(r131)
    r171 = r171 / max(r171)
    r193 = r193 / max(r193)
    r211 = r211 / max(r211)
    r335 = r335 / max(r335)

    resp = [r94, r131, r171, r193, r211, r335]
    linecolors = ['r','g','b', 'c', 'm', 'y']

    fig, axs = plt.subplots()
    # # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    # axs.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

    axs.minorticks_on()
    axs.tick_params(axis='both', which = 'both',
                        direction = 'in',
                        bottom = True, top= True, left= True, right = True)
    axs.tick_params(which = 'major', length = 10)
    axs.tick_params(which = 'minor', length = 5)

    axs.set_xlim(4.5, 7.5)
    axs.set_ylim(0, 1.1)
    axs.set_xlabel(r"$\log_{10}(T/\mathrm{K})$")
    axs.set_ylabel("Normalized AIA Temp. Response")

    for j in range(6):

        respj = resp[j]
        
        # plot responses
        axs.plot(logte, respj, color = linecolors[j], label = allwaves[j])
    
    axs.legend()
    plt.show()

######################################################################
######################################################################


def temp_responses():
    """Plot the temperature response functions for AIA to show where each channel is most sensitive."""

    # load in response functions as 101x1 arrays
    respinfo = sp.io.readsav("aia_response_2011.sav", python_dict = False, verbose = True)
    r94 = respinfo['r94']
    r131 = respinfo['r131']
    r171 = respinfo['r171']
    r193 = respinfo['r193']
    r211 = respinfo['r211']
    r335 = respinfo['r335']
    logte = ((respinfo['r'])['logte'])[0]

    resp = [r94, r131, r171, r193, r211, r335]
    linecolors = ['r','g','b', 'c', 'm', 'y']

    fig, axs = plt.subplots()

    axs.minorticks_on()
    axs.tick_params(axis='both', which = 'both',
                        direction = 'in',
                        bottom = True, top= True, left= True, right = True)
    axs.tick_params(which = 'major', length = 10)
    axs.tick_params(which = 'minor', length = 5)

    axs.set_xlim(4.5, 7.5)
    axs.set_ylim(10**(-28), 10**(-23))
    axs.set_xlabel(r"$\log_{10}(T/\mathrm{K})$")
    axs.set_ylabel(r"AIA Temp. Response $(\mathrm{DN}\;\mathrm{cm}^5 \; \mathrm{px}^{-1})$")

    for j in range(6):
        
        # plot responses
        respj = resp[j]
        axs.semilogy(logte, respj, color = linecolors[j], label = allwaves[j])
    
    axs.legend()
    # plt.style.use("dark_background")
    plt.savefig(fname = "C:\\Users\\Lucien\\Documents\\School\\Research - Thesis\\response.png")#, transparent = True)
