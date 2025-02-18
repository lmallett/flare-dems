# (c) Lucien Mallett 2024
# Produces images from .npz files containing saved variables.

# import matplotlib
# matplotlib.use("Agg")

import matplotlib.colors as colors

import numpy as np
import os
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps.color_tables as aiacolormaps
import astropy.units as u


##############################################
##############################################
# Path to folder containing .npz files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

pathresp = "aia_resp_110809.sav"
path = "D://emcubes_110809//emcubes_np//"
# path = "D://emcubes_no_335_110809//emcubes_no_335_np//"
# dir_emcubes = savedir + "emcubes_no335//"
# xregion = [350,651]
# yregion = [305,746]
xregion = [410,601]
yregion = [410,611]
savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2011//"

##############################################

pathresp = "aia_resp_140910.sav"
path = "D://emcubes_140910//emcubes_np//"
xregion = [250,750]
yregion = [350,600]
savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2014//"

##############################################

dir_emcubes = savedir + "emcubes//"
dir_statuscubes = savedir + "statuscubes//"
dir_satmaps = savedir + "satmaps//"

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

contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

# colors are in comments
color_list = [
    "#000000",  # black, 0 (converges)
    "#984ea3",  # purple, 1 (unbounded)
    "#ffff00",  # yellow, 2 (no soln) 
    "#FF63C2",  # pink, 3 (func did not converge)
    "#576BFF",  # blue,  10 (negative coeff from simplex) 
    "#88DC85",  # green, 11 (nocount; too little data)+
    ]

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
            ax.contour(sat, levels = [0, 1, 3], colors = contourclrs, antialiased = False)#linewidths = 100)
        
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

        fig.savefig(dir_emcubes + file[:-4] + " " + str(i) + ".png", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    dpi = 1,
                    pil_kwargs = {'lossless':True}
                    )        
        
        plt.clf()

        fig = plt.figure(num = 1, clear = True, frameon = False)
        fig.set_size_inches(emcube.shape[1], emcube.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        test_satarr = np.full_like(sat, 0)
        test_satarr[10:25,20:120] = 1
        test_satarr[25:50,40:140] = 2
        test_satarr[50:75,60:160] = 3
        test_satarr[75:100,40:140] = 4
        test_satarr[100:125,20:120] = 5
        test_satarr[125:150,40:140] = 6

        print(test_satarr)
        ax.imshow(test_satarr, cmap = 'gray', interpolation = 'none', origin = 'lower',)
        ax.contour(test_satarr, levels = [0, 1, 3], colors = contourclrs, antialiased = False)#linewidths = 100)
        fig.savefig(dir_emcubes + "test.png", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    dpi = 1,
                    pil_kwargs = {'lossless':True}
                    )
        
        # plt.clf()

    return

# def imgs_emcubes(sat = False, status = False):

#     for i in range(len(filelist)):

#         file = filelist[i]
#         filename = path + file
#         print(file)

#         vars = sp.io.readsav(filename, python_dict = False, verbose = False)
#         emcube = (np.sum(vars.emcube[:,ystart:yend,xstart:xend], axis = 0)**0.25)

#         if sat == TruD:
#             satmap = np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis = 0)
#             satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated
#         if status == TruD:
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
        Saves a series of images into the directory savedir/data_[wave]/.
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
            plt.imsave(savedir + dirs[allwaves[e]] + "a" + str(allwaves[e]) + file[10:-4] + " " + str(i) + ".png",
                data[e,:,:]**0.5,
                cmap = colors[e],
                origin = 'lower')

    plt.close()

def imgs_statuscubes(contour_sat = False):
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

        # 1: The objective function is unbounded.
        # 2: No solution satisfies the given constraints.
        # 3: The routine did not converge.
        statuscube[statuscube == 10] = 4    # SIMPLEX very occasionally returns negative coeffs even though positivity is a constraint. Now these are set to zero and the corresponding status code is 10.
        statuscube[statuscube == 11] = 5    # NOCOUNTS = WHERE(total(image,3) LT 10*eps). Set the status value of nocount pixels to 11.0

        fig = plt.figure(num = 1, clear = True, frameon = False)
        fig.set_size_inches(statuscube.shape[1], statuscube.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        img = ax.imshow(statuscube, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower', interpolation = 'none')

        if contour_sat == True:
            sat = np.sum(np.copy((vars["satmap"])[:,ystart:yend,xstart:xend]), 0)
            sat = sat.astype('int32')
            ax.contour(sat, levels = [0, 2, 4], colors = contourclrs, antialiased = False)#linewidths = 100)

        # satmap = np.sum((vars["satmap"])[:,ystart:yend, xstart:xend], axis = 0)
        # satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated
        # ax.imshow(satmasked, cmap = sat_color, origin = 'lower', interpolation = 'none')

        fig.savefig(dir_statuscubes + file[:-4] + " " + str(i) + ".png", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    dpi = 1,
                    pil_kwargs = {'lossless':True}
                    )
        plt.clf()

    test_arr0 = np.full_like(statuscube, 0)
    test_arr1 = np.full_like(statuscube, 1)
    test_arr2 = np.full_like(statuscube, 2)
    test_arr3 = np.full_like(statuscube, 3)
    test_arr4 = np.full_like(statuscube, 4)
    test_arr5 = np.full_like(statuscube, 5)

    plt.imsave(dir_statuscubes + "test_arr0.png", test_arr0, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower')
    plt.imsave(dir_statuscubes + "test_arr1.png", test_arr1, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower')
    plt.imsave(dir_statuscubes + "test_arr2.png", test_arr2, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower')
    plt.imsave(dir_statuscubes + "test_arr3.png", test_arr3, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower')
    plt.imsave(dir_statuscubes + "test_arr10.png", test_arr4, cmap = color_map, vmin = 0, vmax =5, origin = 'lower')
    plt.imsave(dir_statuscubes + "test_arr11.png", test_arr5, cmap = color_map, vmin = 0, vmax =5, origin = 'lower')

    return


def imgs_satmaps():
    """Produces frames for a movie of the event, which shows how many channels are saturated.
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
        satmap = np.sum((vars["satmap"])[:,ystart:yend, xstart:xend], axis = 0)
        satmap = np.ma.masked_where(satmap == 0, satmap)

        fig = plt.figure(num = 1, clear = True, frameon = False)
        fig.set_size_inches(satmap.shape[1], satmap.shape[0])
        ax = plt.Axes(fig, [0., 0., 1., 1.])
        ax.set_axis_off()
        fig.add_axes(ax)

        cmap = matplotlib.cm.get_cmap("cool").copy()
        cmap.set_bad(color = "black")

        ax.imshow(satmap, cmap = cmap, vmin = 0, vmax = 6, origin = 'lower', interpolation = 'none')

        fig.savefig(dir_satmaps + file[:-4] + ".webp", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    dpi = 1,
                    pil_kwargs = {'lossless':True}
                    )
        
        plt.clf()
    return


def imgs_allchannels():