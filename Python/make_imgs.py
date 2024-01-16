# Lucien Mallett

# This code assumes that, in "savedir/", the following folders exist:
    # - data_94/
    # - data_131/
    # - data_171/
    # - data_193/
    # - data_211/
    # - data_335/
    # - contoured_data_94/
    # - contoured_data_131/
    # - contoured_data_171/
    # - contoured_data_193/
    # - contoured_data_211/
    # - contoured_data_335/
    # - emcubes/
    # - contoured_emcubes/

# and that "path/" contains .sav files, which contain the L1.5 data.

import scipy as sp
import skimage as ski
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import idl_colorbars
import os
import astropy.units as u
import sunpy.visualization.colormaps.color_tables as aiacolormaps

# Path to folder containing .sav files, which contain image data
path = "E:\\2011_08_09\\emcubes\\all\\"

greenwhite = idl_colorbars.getcmap(8)
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
filelist = os.listdir(path)
savedir = "C:/Users/Lucien/Documents/School/Research - Thesis/movies/"

# directories to put individual frames in
dirs = {
    94: "data_94/",
    131: "data_131/",
    171: "data_171/",
    193: "data_193/",
    211: "data_211/",
    335: "data_335/",
    "emcube": "emcubes/"}

xregion = [0,500]
yregion = [0,400]

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

allwaves  = [94,131,171,193,211,335]

######################################################################
######################################################################
######################################################################
######################################################################

def imgs_emcubes():
    """Produces frames for a movie of the event, where the data plotted is the sum of the emission cubes over logT.
    Args:
        None.
    Returns:
        None.
        Saves a series of images into the directory savedir/emcubes/.
    """

    for i in range(1):#len(filelist)):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename, python_dict = False, verbose = False)
        data = (np.sum(vars.emcube, axis = 0)**0.25)

        # plot data
        plt.imshow(data, cmap = greenwhite, origin='lower')
        plt.axis('off')
        plt.imsave(savedir + dirs["emcube"] + file[:-4] + ".png", data, cmap = greenwhite,  origin='lower')

    plt.close()

######################################################################
######################################################################
######################################################################
######################################################################

def imgs_datacubes(wavelength = allwaves):
    """Produces frames for a movie of the event in the specified AIA wavelength (94, 131, 171, 193, 211, or 335).
    Args:
        (Optional) An integer, or list of integers, that is a subset of [94,131,171,193,211,335].
    Returns:
        None.
        Saves a series of images into the directory savedir/data_[wave]/.
    """

    if type(wavelength) != list:
        wavelength = [wavelength]

    imglist = []    # list of datacubes
    for i in range(len(filelist)):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)

        data = np.copy(vars.datacube)
        data[data < 0] = 0
        data = data**0.25

        # save modifed datacubes
        imglist.append(data)

    # maps what wavelength we want to 0-5
    # convenient for accessing datacubes
    indices = [allwaves.index(wave) for wave in wavelength]

    for e in indices:

        colors = aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom) 

        for i in range(len(filelist)):            

            # get data at timestamp, for only the relevant channel
            file = filelist[i]
            data = (imglist[i])[e,:,:]

            plt.imsave(savedir + dirs[allwaves[e]] + file[:-4] + ".png", data, cmap=colors, origin = 'lower')

    plt.close()

######################################################################
######################################################################
######################################################################
######################################################################

def contoured_data(wavelength = allwaves):
    """Produces frames for a movie of the event in the specified AIA wavelength (94, 131, 171, 193, 211, or 335), where saturated pixels are contoured over the data.
    Args:
        (Optional) An integer, or list of integers, that is a subset of [94,131,171,193,211,335].
    Returns:
        None.
        Saves a series of images into the directory savedir/contoured_data_[wave]/ with contours over saturated pixels.
    """
    if type(wavelength) != list:
        wavelength = [wavelength]

    imglist = []
    satimglist = []

    for i in range(len(filelist)):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename)

        data = np.copy(vars.datacube)
        data[data < 0] = 0
        data = data**0.25

        # save modifed datacubes & satmaps
        satimglist.append(vars.satmap)
        imglist.append(data)

    # maps what wavelength we want to 0-5
    # convenient for accessing datacubes
    indices = [allwaves.index(wave) for wave in wavelength]

    for e in indices:

        colors = aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom)

        for i in range(len(filelist)):          

            # get data at timestamp, for only the relevant channel
            file = filelist[i]
            data = (imglist[i])[e,:,:]
            sat = (satimglist[i])[e,:,:]
            print(data.shape)

            fig = plt.figure(frameon = False)
            fig.set_size_inches(data.shape[1], data.shape[0])
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)

            ax.imshow(data, cmap = colors, interpolation='none', origin='lower')
            ax.contour(sat, levels = 1, colors = contourclrs[0], antialiased = False)#linewidths = 100)
            fig.savefig(savedir + "contoured_" + dirs[allwaves[e]] + file[:-4] + ".png", dpi = 1)

            fig.clf()
    plt.close()

            # plt.imshow(data, cmap = colors)
            # plt.contour((satimglist[i])[e,:,:], colors = contourclrs, linewidths = 0.1)
            # plt.axis('off')
            # plt.savefig(savedir + "contoured_" + dirs[allwaves[e]] + file[:-4] + ".png", bbox_inches='tight', pad_inches=0)
            # plt.clf()

######################################################################
######################################################################
######################################################################
######################################################################

def contoured_emcubes():
    """Produces frames for a movie of the event, where the data plotted is the sum of the emission cubes over logT, with saturated pixels contoured on it.
    Args:
        None.
    Returns:
        None.
        Saves a series of images into the directory savedir/contoured_emcubes/ with contours over emission cubes.
    """

    imglist = []
    satimglist = []

    for i in range(len(filelist)):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename)

        sat = np.copy(vars.satmap)
        sat = np.sum(sat, 0)
        sat = sat.astype('int32')

        data = (np.sum(vars.emcube, axis = 0)**0.25)
        print(data.shape)

        # save modifed datacubes & satmaps
        satimglist.append(sat)
        imglist.append(data)
  
    for i in range(len(filelist)):          

        # get data at timestamp, for only the relevant channel
        file = filelist[i]
        data = imglist[i]
        sat  = satimglist[i]

        # this is an awful way to plot the image with contours while retaining pixel sizes, but here we are
        fig = plt.figure(frameon = False)
        fig.set_size_inches(data.shape[1], data.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        ax.imshow(data, cmap = greenwhite, interpolation='none', origin='lower')
        ax.contour(satimglist[i], levels = [1,3,5], colors = contourclrs, antialiased = False)#linewidths = 100)
        fig.savefig(savedir + "contoured_" + dirs["emcube"] + file[:-4] + ".png", dpi = 1)

        fig.clf()
        plt.close()

        # # Attempt involving rounding contours to nearest pixel
        # # https://stackoverflow.com/questions/39642680/create-mask-from-skimage-contour
        # contours = ski.measure.find_contours(sat, level=0.5, fully_connected='low', positive_orientation='low', mask=None)
        # print(contours)
        # clist = np.zeros_like(data, dtype='int')

        # for contour in contours:
        #     clist[np.round(contour[:, 0]).astype('int'), np.round(contour[:, 1]).astype('int')] = 1

        # print(clist)
        # plt.matshow(clist)
        # plt.show()

        # # Without contours
        # x = plt.imshow(data, cmap = greenwhite, interpolation='none')
        # print(type(x))
        # plt.axis('off')
        # plt.savefig(savedir + "contoured_" + dirs["emcube"] + file[:-4] + ".png", format='png', bbox_inches='tight', pad_inches=0)

        # data = np.where(contours==1, A, B)
        # plt.imsave(savedir + "contoured_" + dirs["emcube"] + file[:-4] + ".png", a)

        # fig = plt.figure(frameon = False)
        # fig.set_size_inches(data.shape[1], data.shape[0])
        # ax = plt.Axes(fig, [0., 0., 1., 1.])
        # ax.set_axis_off()
        # fig.add_axes(ax)

        # ax.imshow(data, interpolation='none')
        # fig.savefig(savedir + "contoured_" + dirs["emcube"] + file[:-4] + ".png", dpi = 1)



######################################################################
######################################################################
######################################################################
######################################################################
