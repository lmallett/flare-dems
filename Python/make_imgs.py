# (c) Lucien Mallett 2024
# Produces images from .sav files containing saved variables.

import os

import scipy as sp
import numpy as np

import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import matplotlib.colors as colors
matplotlib.use("Agg")

import astropy.units as u
import sunpy.visualization.colormaps.color_tables as aiacolormaps
import idl_colorbars

##############################################
##############################################
# Path to folder containing .sav files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

# Path to response functions
pathresp = "aia_resp_full.sav"

##############################################

path = "E://emcubes_110809//emcubes_np//"
dir_emcubes = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//emcubes//"
dir_statuscubes = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//statuscubes//"

# path = "E://emcubes_no_335_110809//emcubes_no_335_np//"
# dir_emcubes = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//emcubes_no335//"

# datadir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//"

xregion = [350,651]
yregion = [305,746]

##############################################

path = "E://emcubes_140910//emcubes_np//"
dir_emcubes = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2014//emcubes//"
# datadir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2014//"

xregion = [250,750]
yregion = [350,600]

##############################################

dirs = {
    94: "data_94//",
    131: "data_131//",
    171: "data_171//",
    193: "data_193//",
    211: "data_211//",
    335: "data_335//"
    }

##############################################
##############################################
# Initialization
##############################################
##############################################

greenwhite = idl_colorbars.getcmap(8)
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

color_list = [
    # "#e41a1c", # red
    "#000000",
    # "#4daf4a", # green
    "#61DA5D", # green
    # "#377eb8", # blue
    "#00B8FF", # blue
    # "#984ea3", # purple
    "#FF4FC1", # pink
    # "#ff7f00", # orange
    "#ffff00"    # yellow
    ]

color_list = [
    # "#e41a1c", # red
    "#000000",
    # "#4daf4a", # green
    "#88DC85",
    # "#377eb8", # blue
    # "#33A0F9",
    "#576BFF",
    # "#984ea3", # purple
    # "#FF4FC1", # pink
    "#FF63C2",
    "#ffff00"    # yellow
    ]

# color_list = [
#     '#000000',  # black
#     '#a6cee3',  # l-blue
#     '#1f78b4',  # d-blue
#     '#b2df8a',  # l-green
#     '#33a02c',  # d-green
#     '#eb34b7'   # pink
# ]
# color_list = [
#     "#000000", # black
#     "#40B936", # green
#     "#6BCEFF", # blue
#     "#FF4FC1", # pink
#     "#FFDE3B", # yellow
# ]

color_map = colors.ListedColormap(color_list)
sat_color = colors.ListedColormap("#e41a1c")
# sat_color = colors.ListedColormap("#AF0B00")
# sat_color = colors.ListedColormap("#515477")


allwaves  = [94,131,171,193,211,335]

##############################################
##############################################
# Functions
##############################################
##############################################

def imgs_emcubes(contour_sat = False, color_sat = False, status = False):
    """Produces frames for a movie of the event, where the data plotted is the sum of the emission cubes over logT, with saturated pixels contoured on it.
    Args:
        None.
    Returns:
        None. Saves a series of images into the directory.
    """

    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file
        print(filename)

        vars = np.load(filename)
        emcube = (np.sum( (vars["emcube"])[:,ystart:yend,xstart:xend], axis = 0)**0.25)
        # # could also load in .sav files
        # vars = sp.io.readsav(filename)
        # emcube = (np.sum(vars.emcube[:,ystart:yend,xstart:xend], axis = 0)**0.25)

        fig = plt.figure(num = 1, clear = True, frameon = False)
        fig.set_size_inches(emcube.shape[1], emcube.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        ax.imshow(emcube, cmap = 'gray', interpolation = 'none', origin = 'lower',)

        if contour_sat == True:
            # sat = np.sum(np.copy(vars.satmap[:,ystart:yend,xstart:xend]), 0)
            sat = np.sum(np.copy((vars["satmap"])[:,ystart:yend,xstart:xend]), 0)

            sat = sat.astype('int32')
            ax.contour(sat, levels = [0, 2, 4], colors = contourclrs, antialiased = False)#linewidths = 100)
        
        if color_sat == True:
            # satmap = np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis = 0)
            satmap = np.sum((vars["satmap"])[:,ystart:yend, xstart:xend], axis = 0)

            satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated
            ax.imshow(satmasked, cmap = sat_color, origin = 'lower', interpolation = 'none')

        if status == True: 
            # statuscube = vars.statuscube[ystart:yend, xstart:xend]
            statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]

            statusmasked = np.ma.masked_where(emcube != 0, statuscube)

            ax.imshow(statusmasked, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower', interpolation = 'none')

        fig.savefig(dir_emcubes + file[:-4] + ".webp", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    dpi = 1,
                    pil_kwargs = {'lossless':True}
                    )
        
        plt.clf()
    return

# def imgs_emcubes(sat = False, status = False):

#     for i in range(len(filelist)):

#         file = filelist[i]
#         filename = path + file
#         print(file)

#         vars = sp.io.readsav(filename, python_dict = False, verbose = False)
#         emcube = (np.sum(vars.emcube[:,ystart:yend,xstart:xend], axis = 0)**0.25)

#         if sat == True:
#             satmap = np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis = 0)
#             satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated
#         if status == True:
#             statuscube = vars.statuscube[ystart:yend, xstart:xend]
#             statusmasked = np.ma.masked_where(emcube != 0, statuscube)

#         plt.imsave(dir_emcubes_emcubes + file[:-4] + ".webp", 
#                    data, 
#                    cmap = 'binary_r',  origin='lower',
#                    pil_kwargs={'lossless':True},
#                    format = 'webp')

#     plt.close()

# directories to put individual frames in

def imgs_datacubes(wavelength = allwaves):
    """Produces frames for a movie of the event in the specified AIA wavelength (94, 131, 171, 193, 211, or 335).
    Args:
        (Optional) An integer, or list of integers, that is a subset of [94,131,171,193,211,335].
    Returns:
        None.
        Saves a series of images into the directory datadir/data_[wave]/.
    """

    if type(wavelength) != list:
        print("a")
        wavelength = [wavelength]

    print(wavelength)
    # maps what wavelength we want to indices 0-5, for access purposes
    indices = [allwaves.index(wave) for wave in wavelength]

    colors = [] # list of colormaps
    for e in indices:
        colors.append(aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom))

    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file
        print(filename)

        vars = np.load(filename)
        data = (vars["datacube"])[:, ystart:yend, xstart:xend]
        data[data < 0] = 0

        for e in indices:
            plt.imsave(datadir + dirs[allwaves[e]] + "a" + str(allwaves[e]) + file[10:-4] + ".webp",
                data[e,:,:],
                cmap = colors[e],
                origin = 'lower')

    plt.close()

def imgs_statuscubes():
    """Produces frames for a movie of the event, where the data plotted a color-coded map of the error code
    returned by the DEM solver.
    Args:
        None.
    Returns:
        None. Saves a series of images into the directory.
    """

    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file
        print(filename)

        vars = np.load(filename)

        statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]
        # statusmasked = np.ma.masked_where(emcube != 0, statuscube)

        fig = plt.figure(num = 1, clear = True, frameon = False)
        fig.set_size_inches(statuscube.shape[1], statuscube.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        ax.imshow(statuscube, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower', interpolation = 'none')

        satmap = np.sum((vars["satmap"])[:,ystart:yend, xstart:xend], axis = 0)
        satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated
        ax.imshow(satmasked, cmap = sat_color, origin = 'lower', interpolation = 'none')

        fig.savefig(dir_statuscubes + file[:-4] + ".webp", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    dpi = 1,
                    pil_kwargs = {'lossless':True}
                    )
        
        plt.clf()
    return
