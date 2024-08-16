# IMPORTS

import os
import sys

import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
from matplotlib.backend_bases import MouseButton
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.cm as cm
import idl_colorbars

import matplotlib as mpl
mpl.use('TkAgg')
plt.rcParams['text.usetex'] = True
import matplotlib.gridspec as gridspec

try: 
    mpl.rcParams['keymap.back'].remove('left')
    mpl.rcParams['keymap.forward'].remove('right')
except:
    pass

contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']

##############################################
##############################################
# CHANGE THESE LINES!
##############################################
##############################################

# Path to folder containing .sav files, which contain image data at different timesteps
# path = "E:\\2011_08_09\\emcubes\\all\\"
# path = "E:\\2014_09_10\\emcubes\\all\\"
# path = "E:\\2011_08_09.tar\\2011_08_09\\emcubes\\all\\"


path = "E:\\emcubes_110809\\"
xregion = [350,650]
yregion = [275,775]


# path = "E:\\emcubes_140910\\"
# xregion = [250,750]
# yregion = [350,750]   

# Path to response functions
pathresp = "aia_resp_full.sav"

##############################################
##############################################
##############################################

# load in response functions as 20x1 arrays
respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

logte = (((respinfo['r'])['A94'])[0])['logte'][0]
r94 = (((respinfo['r'])['A94'])[0])['tresp'][0]
r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]

##############################################
##############################################
##############################################

greenwhite = idl_colorbars.getcmap(8)
filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

allwaves  = [94,131,171,193,211,335]

resps = np.array([r94, r131, r171, r193, r211, r335])
linecolors = ['r','g','b', 'c', 'm', 'y']

# Create 2x2 sub plots
gs = gridspec.GridSpec(2, 2)
fig = plt.figure()

ax0 = fig.add_subplot(gs[0, 0]) # row 0, col 0
ax1 = fig.add_subplot(gs[1, 0], sharex = ax0, sharey = ax0) # row 1, col 0
ax2 = fig.add_subplot(gs[:, 1]) # col 1, span all rows

ax0.set_xlim(5.6, 8.0)
ax0.set_ylim(1, 10**6)
ax0.set_yscale("log")
ax0.set_aspect(0.2)
ax0.set_xticks(np.arange(5.6, 8.0, 0.2))
ax0.tick_params(axis='both', which = 'both',
                    direction = 'in',
                    bottom = True, top= True, left= True, right = True)
ax0.tick_params(which = 'major', length = 10)
ax0.tick_params(which = 'minor', length = 5)
ax0.xaxis.set_minor_locator(AutoMinorLocator(2))

fig.set_size_inches(15, 5.5)
fig.tight_layout()
# plt.subplots_adjust(left = 0.07)

# for future blitting optimization
# background1 = fig.canvas.copy_from_bbox(ax2.bbox) 
# background2 = fig.canvas.copy_from_bbox(ax0.bbox) 

image_index = 25

def load_data(i):

    global datacube
    global emcube
    global satmap
    global axs
    global logt
    global sel_ind
    global resp

    # Change what file we're looking at
    file = filelist[i]
    filename = path + file

    vars = sp.io.readsav(filename)
    sel_ind = np.searchsorted(logte, vars.logt)
    logt = logte[sel_ind]
    resp = resps[:,sel_ind]

    datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
    satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
    satmasked = np.ma.masked_where(satmap < 1, satmap) # makes transparent satmap for overplotting
    emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])

    # Plot image
    ax2.clear()
    ax2.imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = 'binary_r', interpolation='none')
    ax2.imshow(satmasked, cmap = 'hsv', origin = 'lower', interpolation = 'none')
    
    # line1 = ax2.axhline(int(datacube.shape[0]/2), color="r", linestyle="--")
    # line2 = ax2.axvline(int(datacube.shape[1]/2), color="r", linestyle="--")

def lociplots(x, y):

    global resp
    global sel_ind

    # causes an ignorable error due to sharey, see https://github.com/matplotlib/matplotlib/issues/9970
    ax0.clear()
    ax1.clear()

    ax0.set_xlim(5.6, 8.0)
    ax0.set_ylim(1, 10**6)

    # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    ax0.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

    # ax2.text(50, 25, f"# saturated channels: {int(satmap[y,x])}", bbox=dict(fill=True, facecolor = 'white', edgecolor='black', linewidth=1))
    ax0.set(title = f"num. saturated channels: {int(satmap[y,x])}")
    ax0.set_xlabel("$log(T)$")
    ax0.set_ylabel("Emission Measure [$10^{26} cm^{-5}$]")

    locis = ( datacube[:,y,x][:,None] / resp ) / (10**26) # produces an array of the 6 functions to be plotted
    loci_avgs = np.average(locis, axis=0).astype(float)
    loci_stdevs= np.std(locis.astype(float), axis=0).astype(float)
    # print(loci_avgs.shape)
    # print(loci_stdevs.shape)
    # print(sp.integrate.simpson(loci_stdevs, x=logt))
    ax0.plot(logt, loci_avgs, color = 'r', linestyle = 'dashed', label = 'loci mean')
    ax1.plot(logt, loci_stdevs, label = 'standard deviation at each temp')
    ax1.grid()
    ax1.set_aspect(0.2)

    #PLOTTING 

    for j in range(6):
        ax0.semilogy(logt, locis[j], color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))        
    
    ax0.legend(loc = 'upper left')
    ax0.xaxis.set_minor_locator(AutoMinorLocator(2))

    fig.show()

def on_move(event):
    if event.inaxes is ax2:

        xcoord = event.xdata
        ycoord = event.ydata
        col = int(xcoord + 0.5)
        row = int(ycoord + 0.5)

        # print(f'data coords row: {row} col: {col}')

        lociplots(col, row)


def on_press(event):
    global image_index

    # print('press', event.key)
    sys.stdout.flush()
    if event.key == "left":
        image_index = (image_index - 1) % len(filelist)
        load_data(image_index)
        fig.canvas.draw()
    elif event.key == "right":
        image_index = (image_index + 1) % len(filelist)
        load_data(image_index)
        fig.canvas.draw()
    print(image_index)

load_data(image_index)

fig.canvas.mpl_connect('key_press_event', on_press)
fig.canvas.mpl_connect('key_press_event', on_move)
fig.canvas.mpl_connect('motion_notify_event', on_move)

plt.show()