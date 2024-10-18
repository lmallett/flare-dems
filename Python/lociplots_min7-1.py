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

import matplotlib.patches as mpatches

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

# path = "E:\\original cadences\\emcubes_110809\\"
# xregion = [350,650]
# yregion = [275,775]
# path = "E:\\emcubes_140910\\"
# xregion = [250,750]
# yregion = [350,750]  

# path = "E:\\emcubes_110809\\emcubes\\"
path = "E:\\emcubes_no_335_110809\\emcubes_no_335\\"
xregion = [400,600]
yregion = [400,600]
# path = "E:\\emcubes_140910\\emcubes\\"
# xregion = [400,600]
# yregion = [400,600]
 

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
outer = gridspec.GridSpec(2, 2, height_ratios=[1,0.1]) # more rows than necessary to control colorbar height
gs0 = gridspec.GridSpecFromSubplotSpec(1, 1, subplot_spec = outer[0])
gs = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec = outer[1], wspace=0.01, width_ratios=[1,1,0.05])
fig = plt.figure()

# ax0 = fig.add_subplot(gs[0, 0]) # row 0, col 0
# ax1 = fig.add_subplot(gs[1, 0], sharex = ax0, sharey = ax0) # row 1, col 0
# ax2 = fig.add_subplot(gs[:, 1]) # col 1, span all rows

# # or dont
# gs = gridspec.GridSpec(1, 2)
# fig = plt.figure()

ax0 = fig.add_subplot(gs0[:, 0]) # row 0, col 0
ax1 = fig.add_subplot(gs[:, 0]) # col 1, span all rows
ax2 = fig.add_subplot(gs[:, 1])
ax2a= fig.add_subplot(gs[:, 2])

ax2.sharey(ax1)

ax1.set_xticks([])
ax1.set_yticks([])
ax2.set_xticks([])
ax2.set_yticks([])

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


# cbar = fig.colorbar(mappable = None, ax = ax2a)
# cax = fig.add_axes([ax2.get_position().x1 + 0.01, ax2.get_position().y0, 0.02, ax2.get_position().height])

fig.set_size_inches(15, 5.5)
# fig.tight_layout()
# plt.subplots_adjust(left = 0.07)

image_index = 30

def load_data(i):

    global datacube
    global emcube
    global satmap
    global logt
    global sel_ind
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

    # confirming which index is logT = 7.1
    # print(logt[15])
    # print(resp[0][15], resp[1][15], resp[-1][15])
    comp_list = np.array([resp[0][15], resp[1][15], resp[-1][15]])

    satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
    emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])
    flatem = np.sum(emcube, axis = 0)
    emmasked = np.ma.masked_where(flatem != 0, flatem) # makes transparent emcube for overplotting

    datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
    data_limited = np.delete(datacube, [2,3,4], axis = 0)   # deletes channels 171, 193, 211


    print("data limited: ", data_limited[:,0,0])
    print("comp list:    ", comp_list)
    respmap = data_limited / comp_list[:, None, None] # gets y/K for logT = 7.1 for 94, 131, 335 channels. [3, x, y]

    respmap = respmap / 1e26
    print("respmap:", respmap[:,0,0])
    print("respmap, full:", (datacube[:,0,0] / resp[:, 15])/ 1e26)

    min_resp = np.min(respmap, axis = 0)                 # finds min of those channels
    arg_min_resp = np.argmin(respmap, axis = 0)
    ratiomap = np.abs(respmap[1,:,:] / min_resp[:,:])     # finds ratio between 131 and the minimum
    ratiomap = np.clip(ratiomap, a_min = 0, a_max= 5)

    print("min_resp:    ", min_resp[0,0])
    print("arg_min_resp:", arg_min_resp[0,0])
    print("ratiomap:    ", ratiomap[0,0])


    # Plot image
    norm = mpl.colors.Normalize(vmin=0, vmax=2)
    cmap0 = mpl.colormaps['binary']
    rgba = [cmap0(0), cmap0(0.5), cmap0(1.0)]
    cmap = mpl.colors.ListedColormap(rgba)

    ax1.clear()
    ax1.imshow(arg_min_resp, origin = 'lower', interpolation='none', cmap = cmap, norm = norm)
    ax1.imshow(emmasked, origin = 'lower', interpolation='none', cmap = 'Reds_r', alpha=0.5)

    patch94 = mpatches.Patch(color= cmap(0), label='Min: 94')
    patch131 = mpatches.Patch(color= cmap(0.5), label='Min: 131')
    patch335 = mpatches.Patch(color= cmap(1.0), label='Min: 335')
    ax1.legend(handles=[patch94, patch131, patch335])

    ax2.clear()
    ratios = ax2.imshow(ratiomap, origin = 'lower', interpolation='none')
    # ax2.imshow(emmasked, origin = 'lower', interpolation='none', cmap = 'Reds_r', alpha=0.5)

    # cbar.remove()
    # ax2a.clear()
    cbar = fig.colorbar(mappable = ratios, cax = ax2a)

    ax1.set_xticks([])
    ax1.set_yticks([])
    ax2.set_xticks([])
    ax2.set_yticks([])

def lociplots(x, y):

    global resp
    global sel_ind
    global ratiomap
    global min_resp
    global arg_min_resp

    # causes an ignorable error due to sharey, see https://github.com/matplotlib/matplotlib/issues/9970
    ax0.clear()

    ax0.set_xlim(5.6, 8.0)
    ax0.set_ylim(1, 10**6)

    # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    ax0.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

    ax0.set(title = f"num. saturated channels: {int(satmap[y,x])}")
    ax0.set_xlabel("$log(T)$")
    ax0.set_ylabel("Emission Measure [$10^{26} cm^{-5}$]")

    locis = ( datacube[:,y,x][:,None] / resp ) / (10**26) # produces an array of the 6 functions to be plotted

    # print(locis[1])
    # resplst = np.array([locis[0][15], locis[1][15], locis[-1][15]])

    # print("Response values:", np.array2string(resplst, formatter = {'float': lambda x: f'{x:10.5f}'}), )
    # print("Minimum:", min_resp[y,x])
    # print("Argminimum:", arg_min_resp[y,x])
    # print("Ratio between 131 and the min:", ratiomap[y,x])
    # print()

    #PLOTTING 

    for j in [0,1,5]:
        ax0.semilogy(logt, locis[j], color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))        
    
    ax0.legend(loc = 'upper left')
    ax0.xaxis.set_minor_locator(AutoMinorLocator(2))
    ax0.vlines(7.1, ymin = 0, ymax = 10**6, colors='r', linestyle = 'dashed')

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