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

import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
from matplotlib.backend_bases import MouseButton
from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)
import matplotlib.cm as cm

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

path = "E:\\emcubes_110809\\emcubes\\"
# path = "E:\\emcubes_no_335_110809\\emcubes_no_335\\"
xregion = [400,600]
yregion = [400,600]
pathresp = "aia_resp_full.sav"

# path = "E:\\emcubes_140910\\emcubes\\"
# xregion = [400,600]
# yregion = [400,600]
# xregion = [250,750]
# yregion = [350,750]  

# Set what image you want to start at
image_index = 45

# And the temperature you want to compare
temp = 6.8
temp = 7.1

# when examining the ratios between the channel and the minimum,
# it can be useful to only look at ratios below a certain value
# specified by this parameter.
clip_value = 5

temp_label = {
    "6.8": 94,
    "7.1": 131
}
temp_idx = {
    "6.8": 0,
    "7.1": 1
}

##############################################
##############################################

filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

# load in response functions as 20x1 arrays
respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

logte= (((respinfo['r'])['A94'])[0])['logte'][0]
r94  = (((respinfo['r'])['A94'])[0])['tresp'][0]
r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]
resps= np.array([r94, r131, r171, r193, r211, r335])

allwaves  = [94,131,171,193,211,335]

linecolors = ['r','g','b', 'c', 'm', 'y']
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
greenwhite = idl_colorbars.getcmap(8)

##############################################
##############################################
# Plotting initialization
##############################################
##############################################

fig = plt.figure()
fig.set_size_inches(15, 5.5)
gs = gridspec.GridSpec(1, 3, width_ratios= [2, 1, 1], figure = fig)

ax0 = fig.add_subplot(gs[0])    # lociplots
ax1 = fig.add_subplot(gs[1])    # left plot
ax2 = fig.add_subplot(gs[2])    # right plot

ax0.set(aspect = 0.2,
        title = "Loci curves, data $y_i$ divided by response function $K_i$",
        xlabel = r"$\mathrm{log(T)}$",
        ylabel = "Emission Measure $[10^{26} \;\mathrm{cm}^{-5}]$",
        xlim = [5.6, 8.0],
        ylim = [1, 1e6],
        xticks = np.arange(5.6, 8.0, 0.2),
        yscale = "log")

ax0.tick_params(axis = 'both', which = 'both', direction = 'in',
                bottom = True, top = True, left = True, right = True)
ax0.tick_params(which = 'major', length = 10)
ax0.tick_params(which = 'minor', length = 5)
ax0.xaxis.set_minor_locator(AutoMinorLocator(2))
ax0.vlines(temp, ymin = 0, ymax = 10**6, colors='r', linestyle = 'dashed')

ax1.set(title = r"Minimum channel at $\mathrm{log}T$" + f"$ = {temp:0.1f}$",
        xticks = [], yticks = [])

norm = mpl.colors.Normalize(vmin = 0, vmax = 2)
cmap0= mpl.colormaps['binary']
rgba = [cmap0(0), cmap0(0.5), cmap0(1.0)]
cmap = mpl.colors.ListedColormap(rgba)
patch94  = mpatches.Patch(facecolor = cmap(0.0), edgecolor = 'k', label='Min: 94')
patch131 = mpatches.Patch(facecolor = cmap(0.5), edgecolor = 'k', label='Min: 131')
patch335 = mpatches.Patch(facecolor = cmap(1.0), edgecolor = 'k', label='Min: 335')
ax1.legend(handles=[patch94, patch131, patch335],
        #    bbox_to_anchor=(1, 0)
            )

ax2.sharey(ax1)
ax2.set(title = f"Ratio of $y_{{{temp_label[str(temp)]}}} / K_{{{temp_label[str(temp)]}}}$ to minimum $y_i/K_i$",
        xticks = [], yticks = [])

##############################################
##############################################
# Plotting interaction
##############################################
##############################################

def load_data(i):

    global datacube
    global emcube
    global satmap

    global logt
    global sel_ind
    global temp_index
    global resp
    global cbar

    global ratiomap
    global min_resp
    global arg_min_resp

    # Change what file we're looking at
    file = filelist[i]
    filename = path + file

    vars = sp.io.readsav(filename)
    sel_ind = np.searchsorted(logte, vars.logt)
    logt = logte[sel_ind]
    resp = resps[:,sel_ind]

    # confirming which index is logT = temp
    temp_index = (np.where(logt == temp))[0][0]     # np returns tuple so we need to index
    comp_list = np.array([resp[0][temp_index], resp[1][temp_index], resp[-1][temp_index]])
    print(logt[temp_index])
    print(comp_list)

    satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
    emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])
    flatem = np.sum(emcube, axis = 0)
    emmasked = np.ma.masked_where(flatem != 0, flatem)      # make transparent emcube for overplotting

    datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
    data_limited = np.delete(datacube, [2,3,4], axis = 0)   # deletes channels 171, 193, 211

    respmap = data_limited / comp_list[:, None, None]       # gets y/K for logT = 7.1 for 94, 131, 335 channels. [3, x, y]
    respmap = respmap / 1e26

    min_resp = np.min(respmap, axis = 0)                # finds min of those channels
    arg_min_resp = np.argmin(respmap, axis = 0)
    
    print(temp_idx[str(temp)])
    ratiomap = np.abs(respmap[temp_idx[str(temp)],:,:] / min_resp[:,:])   # finds ratio between 131 and the minimum
    ratiomap = np.clip(ratiomap, a_min = 0, a_max= 5)

    # Plot image
    ax1.imshow(arg_min_resp, origin = 'lower', interpolation = 'none', cmap = cmap, norm = norm)
    ax1.imshow(emmasked, origin = 'lower', interpolation = 'none', cmap = 'Reds_r', alpha= 0.5)

    ratios = ax2.imshow(ratiomap, origin = 'lower', interpolation ='none', cmap = 'viridis', vmin = 0, vmax = clip_value)
    # ax2.imshow(emmasked, origin = 'lower', interpolation='none', cmap = 'Reds_r', alpha=0.5)

    cax = fig.add_axes([ax2.get_position().x1 + 0.01, ax2.get_position().y0, 0.01, ax2.get_position().height])
    cbar = plt.colorbar(mappable = ratios, cax = cax,
                ticks = np.arange(0, clip_value+1)
                )

def lociplots(x, y):

    global resp
    global sel_ind
    global temp_index
    global ratiomap

    global min_resp
    global arg_min_resp

    # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    locis = ( datacube[:,y,x][:,None] / resp ) / (10**26) # produces an array of the 6 functions to be plotted

    # Plotting
    for line in list(ax0.lines):
        line.remove()
    for txt in list(ax0.texts):
        txt.remove()

    ax0.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')
    ax0.text(0, -0.2, f"num. saturated channels: {int(satmap[y,x])}", fontsize = 'medium', transform=ax0.transAxes)

    for j in [0,1,5]:
        ax0.semilogy(logt, locis[j], color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))        
    ax0.legend(loc = 'upper left')

    fig.show()

def on_move(event):
    if event.inaxes is ax1 or event.inaxes is ax2:

        xcoord = event.xdata
        ycoord = event.ydata
        col = int(xcoord + 0.5)
        row = int(ycoord + 0.5)

        # print(f'data coords row: {row} col: {col}')
        lociplots(col, row)

def on_press(event):
    global image_index

    global cbar
    cbar.remove()

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
    if event.inaxes is ax1 or event.inaxes is ax2:
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