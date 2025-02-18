# (c) Lucien Mallett 2024
# This program creates an interactive plot.
# On the right panel is the summed emission cube.
# On the left panel is the loci curves (data divided by response function, for each channel).
# When you mouse over the right panel, the loci plots for that particular pixel
# are shown on the left.
# Use left/right arrow key to go to the prev/next emcube files.

import os
import sys

import numpy as np
import scipy as sp
import idl_colorbars

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.ticker import AutoMinorLocator
import matplotlib.colors as colors

from textwrap import wrap

mpl.use('TkAgg')
plt.rcParams['text.usetex'] = True

try: 
    mpl.rcParams['keymap.back'].remove('left')
    mpl.rcParams['keymap.forward'].remove('right')
except:
    pass

##############################################
##############################################
# Path to folder containing .sav files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

# Path to response functions
pathresp = "aia_resp_110809.sav"

path = "D:\\emcubes_110809\\emcubes\\"
# path = "D:\\emcubes_no_335_110809\\emcubes_no_335\\"
# xregion = [350,650]
# yregion = [275,775]
xregion = [410,601]
yregion = [410,611]

pathresp = "aia_resp_140910.sav"
path = "D:\\emcubes_140910\\emcubes\\"
# path = "D:\\emcubes_140910\\emcubes_no_335\\"
xregion = [250,750]
yregion = [360,600]   

xregion = [375,575]
yregion = [400,480]

# Set what image you want to start at
image_index = 45

##############################################
##############################################
# Initialization
##############################################
##############################################

filelist = os.listdir(path)

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

allwaves = [94,131,171,193,211,335]

resps = np.array([r94, r131, r171, r193, r211, r335])

linecolors = ['r','g','b', 'c', 'm', 'y']
# contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
greenwhite = idl_colorbars.getcmap(8)

# color_map1 = ['#000000', '#a6cee3', '#1f78b4', '#b2df8a', '#33a02c', '#eb34b7']
# color_list = ['#fdbf6f', '#a6cee3', '#1f78b4', '#b2df8a', '#eb34b7']
# color_map = colors.ListedColormap(color_map1) 

color_list = [
    "#e41a1c", # red
    "#4daf4a", # green
    "#377eb8", # blue
    "#984ea3", # purple
    "#ff7f00", # orange
    ]
color_map = colors.ListedColormap(color_list)
sat_color = colors.ListedColormap(color_list[0])

# clean handling of multiline custom legend
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
# Plotting initialization
##############################################
##############################################

fig = plt.figure()
fig.set_size_inches(12, 4.5)
gs = gridspec.GridSpec(1, 2, wspace= 0.1)
plt.subplots_adjust(left = 0.05, right=0.85)

ax0 = fig.add_subplot(gs[:, 0]) # locitplots
ax1 = fig.add_subplot(gs[:, 1]) # image plot

ax0.set(aspect = 0.2,
        title = "Loci curves, $y_i / K_i$",
        xlabel = r"$\mathrm{log(T)}$",
        ylabel = "Emission Measure $[10^{26} \;\mathrm{cm}^{-5}]$",
        xlim = [5.6, 8.0],
        ylim = [1, 1e6],
        xticks = np.arange(5.6, 8.0, 0.2),
        yscale = "log")

ax0.tick_params(axis='both', which = 'both',
                    direction = 'in',
                    bottom = True, top= True, left= True, right = True)
ax0.tick_params(which = 'major', length = 10)
ax0.tick_params(which = 'minor', length = 5)
ax0.xaxis.set_minor_locator(AutoMinorLocator(2))

##############################################
##############################################
# Plotting interaction
##############################################
##############################################

def load_data(i):
    """Produces the image for the current file that is loaded -- sets up right-hand panel."""
    
    global datacube
    global emcube
    global statuscube
    global satmap

    global logt
    global resp
    global sel_ind

    # Change what file we're looking at
    file = filelist[i]
    filename = path + file

    vars = sp.io.readsav(filename)
    sel_ind = np.searchsorted(logte, vars.logt)
    logt = logte[sel_ind]
    resp = resps[:,sel_ind]

    # IDL is col-major, so the data is unfortunately stored as [channel, y, x]
    datacube = vars.datacube[:,ystart:yend, xstart:xend]

    satmap = np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis = 0)
    satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated

    emcube = vars.emcube[:,ystart:yend, xstart:xend]
    emcubesummed = np.sum(emcube, axis = 0)**0.25

    # makes array describing why each DEM failed
    statuscube = vars.statuscube[ystart:yend, xstart:xend]
    statusmasked = np.ma.masked_where(emcubesummed != 0, statuscube)

    # Plot image
    ax1.clear()
    ax1.imshow(emcubesummed, origin = 'lower', cmap = 'binary_r', interpolation='none')
    ax1.text(0, -0.2, f"{file[11:-4]}", fontsize = 'medium', transform=ax1.transAxes)

    # choose cmap = 'hsv' or other map to see how saturated certain areas are
    ax1.imshow(statusmasked, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower', interpolation = 'none')
    ax1.imshow(satmasked, cmap = sat_color, origin = 'lower', interpolation = 'none')
    ax1.set_title("Error Codes")

    ax1.legend(handles = patches,
               bbox_to_anchor=(1.05, 1),
               loc = 2,
               borderaxespad = 0,
               fontsize = 'small',
               labelspacing = 1)

def lociplots(x, y):
    """Handles left-hand plot based on x, y position of the mouse on the right-hand side of the plot.
    For each channel, divide the data by its corresponding response function.
    """

    global resp
    global sel_ind
    global statuscube

    locis = ( datacube[:,y,x][:,None] / resp ) / (10**26) # produces an array of the 6 functions to be plotted

    # sometimes causes an ignorable error due to sharey
    # see https://github.com/matplotlib/matplotlib/issues/9970
    # ax1.clear()

    # Plotting 
    for line in list(ax0.lines):
        line.remove()
    for txt in list(ax0.texts):
        txt.remove()

    ax0.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')
    ax0.text(0, -0.2, f"num. saturated channels: {int(satmap[y,x])}", fontsize = 'medium', transform=ax0.transAxes)

    for j in range(6):
        ax0.semilogy(logt, locis[j], color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))            
    ax0.legend(loc = 'upper left')

    fig.show()

def on_move(event):
    """Updates loci plots with mouse movement."""

    if event.inaxes is ax1:

        xcoord = event.xdata
        ycoord = event.ydata
        col = int(xcoord + 0.5)
        row = int(ycoord + 0.5)

        # print(f'data coords row: {row} col: {col}')
        lociplots(col, row)

def on_press(event):
    """Handles flipping between images."""

    global image_index

    sys.stdout.flush()
    if event.key == "left":
        image_index = (image_index - 1) % len(filelist)
        load_data(image_index)
        fig.canvas.draw()
    elif event.key == "right":
        image_index = (image_index + 1) % len(filelist)
        load_data(image_index)
        fig.canvas.draw()
    print(filelist[image_index])

    # also updates the lociplots if necessary
    if event.inaxes is ax1:
        xcoord = event.xdata
        ycoord = event.ydata
        col = int(xcoord + 0.5)
        row = int(ycoord + 0.5)
        lociplots(col, row)

load_data(image_index)

fig.canvas.mpl_connect('key_press_event', on_press)
fig.canvas.mpl_connect('key_press_event', on_move)
fig.canvas.mpl_connect('motion_notify_event', on_move)

plt.show()