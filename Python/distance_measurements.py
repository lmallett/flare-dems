# (c) Lucien Mallett 2024

import scipy as sp
import numpy as np
import os

import matplotlib
# matplotlib.use("Agg")
matplotlib.use('TkAgg')

from matplotlib import rc
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import idl_colorbars
import astropy.units as u
import sunpy.visualization.colormaps.color_tables as aiacolormaps

#####################################################
# Path to folder containing .sav files which contain image data

# OLD LOCATIONS
# path = "E:\\2011_08_09\\emcubes\\all\\"
# path = "E:\\2014_09_10\\emcubes\\all\\"
# path = "E:\\emcubes_110809\\uncropped\\"
# path = "E:\\emcubes_140910\\uncropped\\"

path = "E:\\emcubes_110809\\"
xregion = [350,650]
yregion = [275,775]

# path = "E:\\emcubes_140910\\"
# xregion = [250,750]
# yregion = [350,750]   

# Path to response functions 
pathresp = "aia_resp_full.sav"

# Main save folder
savedir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//"


x = 144
y = 214
# x = 250
# y = 144

#####################################################

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

allwaves  = [94,131,171,193,211,335]
respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

logte = (((respinfo['r'])['A94'])[0])['logte'][0]
r94 = (((respinfo['r'])['A94'])[0])['tresp'][0]
r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]

greenwhite = idl_colorbars.getcmap(8)
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
filelist = os.listdir(path)

#####################################################

def r(wavelength, x):
    respcall = "r" + str(wavelength) + "[x]"
    print(exec(respcall))
    return exec(respcall)

#####################################################

resp = np.array([r94, r131, r171, r193, r211, r335])

linecolors = ['r','g','b', 'c', 'm', 'y']
fig, axs = plt.subplots(1,2)
fig.set_size_inches(11, 5.5)
fig.tight_layout()
plt.subplots_adjust(left = 0.07)

# FOR EACH IMAGE
for i in range(1): #len(filelist)):

    # Create 2x2 sub plots
    gs = gridspec.GridSpec(2, 2)
    fig = plt.figure()
    ax0 = fig.add_subplot(gs[0, 0]) # row 0, col 0
    ax1 = fig.add_subplot(gs[1, 0], sharex = ax0, sharey = ax0) # row 1, col 0
    ax2 = fig.add_subplot(gs[:, 1]) # col 1, span all rows

    file = filelist[i]
    filename = path + file
    print(file)

    vars = sp.io.readsav(filename)
    datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
    emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])

    # only select indices at temperature bins
    sel_ind = np.searchsorted(logte, vars.logt)
    logt = logte[sel_ind]

    resp = resp[:,sel_ind]
    locis = ( datacube[:,y,x][:,None] / resp ) / (10**26) # produces an array of the 6 functions to be plotted


    loci_avgs = np.average(locis, axis=0).astype(float)
    loci_stdevs= np.std(locis.astype(float), axis=0).astype(float)
    print(loci_avgs.shape)
    print(loci_stdevs.shape)
    print(sp.integrate.simpson(loci_stdevs, x=logt))

    #PLOTTING 

    for j in range(6):
        ax0.semilogy(logt, locis[j], color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))        
    
    ax0.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')
    ax0.tick_params(axis='both', which = 'both',
                        direction = 'in',
                        bottom = True, top= True, left= True, right = True)
    ax0.tick_params(which = 'major', length = 10)
    ax0.tick_params(which = 'minor', length = 5)
    ax0.set_xlim(5.6, 8.0)
    ax0.set_ylim(1, 10**6)
    ax0.set_xlabel("$log(T)$")
    ax0.set_ylabel("Emission Measure [$10^{26} \mathrm{cm}^{-5}$]")
    ax0.set_yscale("log")
    ax0.set_aspect(0.2)
    ax0.set_xticks(np.arange(5.6, 8.0, 0.2))
    ax0.legend(loc = 'upper left')

    ax0.plot(logt, loci_avgs, color = 'r', linestyle = 'dashed')

    ax1.plot(logt, loci_stdevs)
    ax1.grid()
    ax1.set_aspect(0.2)

    ax2.imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = 'binary_r')
    ax2.axhline(y, color="r", linestyle="--")
    ax2.axvline(x, color="r", linestyle="--")
    ax2.axis('off')
    fig.show()
