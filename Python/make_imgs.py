import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import idl_colorbars
import os
import astropy.units as u
import sunpy.visualization.colormaps.color_tables as aiacolormaps
import types

# Path to folder containing .sav files, which contain image data
path = "E:\\2011_08_09\\emcubes\\all\\"

greenwhite = idl_colorbars.getcmap(8)
filelist = os.listdir(path)
savedir = "C:/Users/Lucien/Documents/School/Research - Thesis/movies/"


xregion = [0,500]
yregion = [0,400]

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1


allwaves  = [94,131,171,193,211,335]

def imgs_emcubes():
    """Produces movies of the event, where data is the sum of the emission cubes.

    Args:
        None.
        
    Returns:
        None.
        Produces a video "datacubes[wavelength].mp4" in the same file folder that is a video of the event,
        in each wavelength specified.
    """

    fig, ax = plt.subplots()

    imglist = []

    for i in range(len(filelist)):
        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename, python_dict = False, verbose = False)
        print(vars.keys())
        print(np.info(vars.datacube))

        data = (np.sum(vars.emcube, axis = 0)**0.25)

        frame = ax.imshow(data, cmap = greenwhite)     # can also use matshow
        text = ax.text(0.5,1.05,file, 
                        size = plt.rcParams["axes.titlesize"],
                        ha = "center",
                        transform = ax.transAxes,
                        backgroundcolor = 'white')
        
        if i == 0:                          # show an initial one first
            ax.imshow(data, cmap = greenwhite, interpolation = "nearest")
            ax.text(0.5,1.05,file, 
                        size = plt.rcParams["axes.titlesize"],
                        ha = "center",
                        transform = ax.transAxes,
                        backgroundcolor = 'white')
            
        imglist.append([frame, text]) 

    ani = animation.ArtistAnimation(fig, imglist, interval = 50, blit = False, repeat_delay = 500)
    ani.save("emcubes.mp4")
    plt.show()


    # Another attempt at using MPL to create images

    # def updatefig(i):
    #     file = filelist[i]
    #     filename = path + file
    #     vars = sp.io.readsav(path + file, python_dict=False, verbose=False)
    #     a = np.sum(vars.emcube, axis = 0)**0.25
    #     frame = ax.imshow(a, cmap = greenwhite)
    #     text = ax.text(0.5,1.05,file, 
    #                     size=plt.rcParams["axes.titlesize"],
    #                     ha="center",
    #                     transform=ax.transAxes,
    #                     backgroundcolor = 'white')
    #     return frame, text
    # ani = animation.FuncAnimation(fig, updatefig, interval = 50, blit = True)



def imgs_datacubes(wavelength = allwaves):
    """Produces movies of the event in the specified AIA wavelength (94, 131, 171, 193, 211, or 335).

    Args:
        (Optional) An integer, or list of integers, that is a subset of [94,131,171,193,211,335].

    Returns:
        None.
        Produces a video "datacubes[wavelength].mp4" in the same file folder that is a video of the event,
        in each wavelength specified.
    """

    if type(wavelength) != list:
        wavelength = [wavelength]

    fig, ax = plt.subplots()

    imglist = []
    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)

        data = np.copy(vars.datacube)
        data[data < 0] = 0
        data = data**0.25

        imglist.append(data)

    indices = [allwaves.index(wave) for wave in wavelength] # maps what wavelength we want to 0-5

    for e in indices:
        colors = aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom)  # set colormap
        movie = []
        for i in range(len(filelist)):            

            file = filelist[i]
            data = (imglist[i])[e,:,:]          # get data at timestamp, for only the relevant channel

            frame = ax.imshow(data, cmap = colors)     # can also use matshow, default = greenwhite
            text = ax.text(0.5,1.05,file, 
                            size = plt.rcParams["axes.titlesize"],
                            ha = "center",
                            transform = ax.transAxes,
                            backgroundcolor = 'white')
            plt.axis('off')
            
            if i == 0:                          # show an initial one first
                ax.imshow(data, cmap = colors, interpolation = "nearest")
                ax.text(0.5,1.05,file, 
                            size = plt.rcParams["axes.titlesize"],
                            ha = "center",
                            transform = ax.transAxes,
                            backgroundcolor = 'white')
                plt.axis('off')
            movie.append([frame, text]) 

        ani = animation.ArtistAnimation(fig, movie, interval = 50, blit = False, repeat_delay = 500)
        ani.save(savedir + "/datacubes_" + str(allwaves[e]) + ".mp4")
        # plt.show()
        # plt.close()


def contoured_data(wavelength = allwaves):

    if type(wavelength) != list:
        wavelength = [wavelength]

    fig, ax = plt.subplots()

    imglist = []
    satimglist = []

    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file

        print(file)
        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)

        data = np.copy(vars.datacube)
        data[data < 0] = 0
        data = data**0.25

        satimglist.append(vars.satmap)
        imglist.append(data)

    indices = [allwaves.index(wave) for wave in wavelength] # maps what wavelength we want to 0-5

    for e in indices:
        colors = aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom)  # set colormap
        movie = []
        for i in range(len(filelist)):            

            file = filelist[i]
            data = (imglist[i])[e,:,:]                  # get data at timestamp, for only the relevant channel

            if i == 0:                          # show an initial one first
                ax.imshow(data, cmap = colors, interpolation = "nearest")
                ax.text(0.5,1.05,file, 
                            size = plt.rcParams["axes.titlesize"],
                            ha = "center",
                            transform = ax.transAxes,
                            backgroundcolor = 'white')
                ax.contour((satimglist[i])[e,:,:])

            frame = ax.imshow(data, cmap = colors)      # can also use matshow, default = greenwhite
            text = ax.text(0.5,1.05,file, 
                            size = plt.rcParams["axes.titlesize"],
                            ha = "center",
                            transform = ax.transAxes,
                            backgroundcolor = 'white')
            contour = plt.contour((satimglist[i])[e,:,:], colors = 'black', linewidths = 0.1)
        
            movie.append(contour.collections + [frame, text]) 


        ani = animation.ArtistAnimation(fig, movie, interval = 50, blit = False, repeat_delay = 500)
        ani.save(savedir + "datacubes_contoured_" + str(allwaves[e]) + ".mp4")
        plt.show()