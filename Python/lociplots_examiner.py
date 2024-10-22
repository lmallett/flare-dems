# IMPORTS

import os
import sys

import matplotlib.pyplot as plt

import numpy as np
import scipy as sp
import idl_colorbars

import matplotlib as mpl
from matplotlib.ticker import (AutoMinorLocator)
import matplotlib.colors as colors
import matplotlib.patches as mpatches
from textwrap import wrap

mpl.use('TkAgg')
plt.rcParams['text.usetex'] = True
import matplotlib.gridspec as gridspec

try: 
    mpl.rcParams['keymap.back'].remove('left')
    mpl.rcParams['keymap.forward'].remove('right')
except:
    pass

##############################################
##############################################
# Path to folder containing .sav files, which contain image data at different timesteps
##############################################
##############################################

path = "E:\\emcubes_110809\\emcubes\\"
# path = "E:\\emcubes_no_335_110809\\emcubes_no_335\\"
xregion = [350,650]
yregion = [275,775]

# path = "E:\\emcubes_140910\\emcubes\\"
# path = "E:\\emcubes_140910\\emcubes_no_335\\"
# xregion = [250,750]
# yregion = [350,750]   

##############################################
##############################################

# Path to response functions
pathresp = "aia_resp_full.sav"

# load in response functions as 20x1 arrays
respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

logte= (((respinfo['r'])['A94'])[0])['logte'][0]
r94  = (((respinfo['r'])['A94'])[0])['tresp'][0]
r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]

##############################################
##############################################

linecolors = ['r','g','b', 'c', 'm', 'y']

# contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
# greenwhite = idl_colorbars.getcmap(8)

# color_map1 = ['#000000', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#eb34b7']
# color_map = colors.ListedColormap(color_map1) 

# color_list = ['#fdbf6f', '#a6cee3', '#1f78b4', '#b2df8a', '#eb34b7']

color_list = [
    "#e41a1c", # red
    "#4daf4a", # green
    "#377eb8", # blue
    "#984ea3", # purple
    "#ff7f00", # orange
    ]
color_map = colors.ListedColormap(color_list)
sat_color = colors.ListedColormap(color_list[0])

labels = [
        "Data is saturated.",
        "1: The objective function is unbounded.",
        "2: No solution satisfies the given constraints.",
        "3: The routine did not converge.",
        "10: Solver sometimes returns negative coeffs in spite of constraints."
        ]

labels = ['\n'.join(wrap(v,20)) for v in labels]
patches =[mpatches.Patch(color = color_list[i], label = labels[i]) for i in range(len(color_list))]

##############################################
##############################################

filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

allwaves  = [94,131,171,193,211,335]

resps = np.array([r94, r131, r171, r193, r211, r335])

##############################################
##############################################

gs = gridspec.GridSpec(1, 2, wspace= 0.01)
fig = plt.figure()

ax0 = fig.add_subplot(gs[0, 0]) # row 0, col 0
ax1 = fig.add_subplot(gs[:, 1]) # col 1, span all rows

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
# fig.tight_layout()
# plt.subplots_adjust(left = 0.07)

image_index = 35

def load_data(i):

    global datacube
    global emcube
    global satmap
    global logt
    global sel_ind
    global resp
    global statuscube

    # Change what file we're looking at
    file = filelist[i]
    filename = path + file

    vars = sp.io.readsav(filename)
    sel_ind = np.searchsorted(logte, vars.logt)
    logt = logte[sel_ind]
    resp = resps[:,sel_ind]

    datacube = vars.datacube[:,ystart:yend, xstart:xend]
    satmap = np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0)
    satmasked = np.ma.masked_where(satmap < 1, satmap) # makes transparent satmap for overplotting
    emcube = vars.emcube[:,ystart:yend, xstart:xend]
    emcubesummed = np.sum(emcube, axis = 0)**0.25

    statuscube = vars.statuscube[ystart:yend, xstart:xend]
    statusmasked = np.ma.masked_where(emcubesummed != 0, statuscube)

    # Plot image
    ax1.clear()
    ax1.imshow(emcubesummed, origin = 'lower', cmap = 'binary_r', interpolation='none')

    # choose cmap = 'hsv' or other map to see how saturated certain areas are
    ax1.imshow(statusmasked, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower', interpolation = 'none')
    ax1.imshow(satmasked, cmap = sat_color, origin = 'lower', interpolation = 'none')

    ax1.legend(handles = patches,
               bbox_to_anchor=(1.05, 1),
               loc = 2,
               borderaxespad = 0,
               fontsize = 'small',
               labelspacing = 1)

def lociplots(x, y):

    global resp
    global sel_ind
    global statuscube

    # sometimes causes an ignorable error due to sharey, see https://github.com/matplotlib/matplotlib/issues/9970
    ax0.clear()

    ax0.set_xlim(5.6, 8.0)
    ax0.set_ylim(1, 10**6)

    # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    ax0.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

    ax0.set(title = f"num. saturated channels: {int(satmap[y,x])} \n status code: {int(statuscube[y,x])}")

    ax0.set_xlabel("$log(T)$")
    ax0.set_ylabel("Emission Measure [$10^{26} cm^{-5}$]")

    locis = ( datacube[:,y,x][:,None] / resp ) / (10**26) # produces an array of the 6 functions to be plotted

    #PLOTTING 

    for j in range(6):
        ax0.semilogy(logt, locis[j], color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))        
    
    ax0.legend(loc = 'upper left')
    ax0.xaxis.set_minor_locator(AutoMinorLocator(2))

    fig.show()

def on_move(event):
    if event.inaxes is ax1:

        xcoord = event.xdata
        ycoord = event.ydata
        col = int(xcoord + 0.5)
        row = int(ycoord + 0.5)

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