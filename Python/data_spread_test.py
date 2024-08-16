# (c) Lucien Mallett 2024

import scipy as sp
import numpy as np
import os
import sys

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

filelist = os.listdir(path)

#####################################################

linecolors = ['r','g','b', 'c', 'm', 'y']
fig, axs = plt.subplots(1,2)
fig.tight_layout()
plt.subplots_adjust(left = 0.07)

image_index = 25

def load_data(i):

    plt.cla()
    file = filelist[i]
    filename = path + file
    print(file)

    vars = sp.io.readsav(filename)
    datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
    emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])
    emcubesum = np.sum(emcube, axis = 0)**0.25

    data_stdevs= np.std(datacube.astype(float), axis=0).astype(float)

    print(data_stdevs.shape)

    axs[0].imshow(emcubesum, origin = 'lower', cmap = 'binary_r')
    axs[1].imshow(data_stdevs, origin = 'lower', cmap = 'binary_r')
    axs[1].contour(emcubesum, [0])
    fig.show()

def on_press(event):
    global image_index

    print('press', event.key)
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