# IMPORTS

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from matplotlib.backend_bases import MouseButton
import matplotlib.cm as cm
import idl_colorbars

import matplotlib as mpl
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
path = "E:\\2011_08_09\\emcubes\\all\\"
# path = "E:\\2014_09_10\\emcubes\\all\\"
# path = "E:\\2011_08_09.tar\\2011_08_09\\emcubes\\all\\"

# Path to response functions
pathresp = "aia_response_2011.sav"

##############################################
##############################################
##############################################

greenwhite = idl_colorbars.getcmap(8)
filelist = os.listdir(path)

xregion = [0,500]
yregion = [0,400]

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

allwaves  = [94,131,171,193,211,335]

# load in response functions as 20x1 arrays
respinfo = sp.io.readsav(pathresp, python_dict = False, verbose = False)
r94 = respinfo['r94']
r131 = respinfo['r131']
r171 = respinfo['r171']
r193 = respinfo['r193']
r211 = respinfo['r211']
r335 = respinfo['r335']
logt = respinfo['logt']
logte = ((respinfo['r'])['logte'])[0]

# only select data at the same logT values (response functions are defined at more logT values than DEM solutions)
sel_ind = respinfo['sel_ind']
resp = [r94, r131, r171, r193, r211, r335]
linecolors = ['r','g','b', 'c', 'm', 'y']

# Set up window
fig, axs = plt.subplots(1,2)
fig.set_size_inches(11, 5.5)
fig.tight_layout()
plt.subplots_adjust(left = 0.07)

axs[1].clear()
axs[1].set_xlim(5.6, 7.5)
axs[1].set_ylim(1, 10**6)
axs[1].set_xlabel("$log(T)$")
axs[1].set_ylabel("Emission Measure [$10^{26} cm^{-5}$]")
axs[1].set_yscale("log")
axs[1].set_aspect(0.2)
axs[1].set_xticks(np.arange(5.6, 7.5, 0.2))
axs[1].tick_params(axis='both', which = 'both',
                    direction = 'in',
                    bottom = True, top= True, left= True, right = True)
axs[1].tick_params(which = 'major', length = 10)
axs[1].tick_params(which = 'minor', length = 5)

# for future blitting optimization
# background1 = fig.canvas.copy_from_bbox(axs[0].bbox) 
# background2 = fig.canvas.copy_from_bbox(axs[1].bbox) 

image_index = 0

def load_data(i):

    global datacube
    global emcube
    global satmap
    global axs

    # Change what file we're looking at
    file = filelist[i]
    filename = path + file

    vars = sp.io.readsav(filename)
    print(vars.logt)
    datacube = np.copy(vars.datacube)
    satmap = np.copy(np.sum(vars.satmap, axis=0))
    satmasked = np.ma.masked_where(satmap < 1, satmap) # makes transparent satmap for overplotting
    emcube = np.copy(vars.emcube)

    # Plot image
    axs[0].clear()
    axs[0].imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = 'binary_r', interpolation='none')
    axs[0].imshow(satmasked, cmap = 'hsv', origin = 'lower', interpolation = 'none')


def lociplots(x, y):

    axs[1].clear()
    axs[1].set_xlim(5.6, 7.5)
    axs[1].set_ylim(1, 10**6)

    # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    axs[1].step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

    # axs[0].text(50, 25, f"# saturated channels: {int(satmap[y,x])}", bbox=dict(fill=True, facecolor = 'white', edgecolor='black', linewidth=1))
    axs[1].set(title = f"# saturated channels: {int(satmap[y,x])}")

    for j in range(6):

        respj = resp[j]
        respj = respj[sel_ind] # 20x1 array

        # plot responses
        # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
        datarray = np.full_like(respj, datacube[j,y,x])
        datarray2 = datarray / respj
        datarray2 = ( datarray2/(10**(26)) ) # 20x1 array

        axs[1].semilogy(logt, datarray2, color = linecolors[j], label = allwaves[j])

    axs[1].legend(loc = 'upper left')
    fig.show()

def on_move(event):
    if event.inaxes:

        xcoord = event.xdata
        ycoord = event.ydata
        col = int(xcoord + 0.5)
        row = int(ycoord + 0.5)

        print(f'data coords row: {row} col: {col}')

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