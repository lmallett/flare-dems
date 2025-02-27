# (c) Lucien Mallett 2024
# Produces images from .npz files containing saved variables.

import matplotlib
matplotlib.use("Agg")

import matplotlib.colors as colors
import matplotlib as mpl

import numpy as np
import os
import matplotlib.pyplot as plt
import sunpy.visualization.colormaps.color_tables as aiacolormaps
import astropy.units as u


plt.rc('text',usetex = True)
font = {'family':'serif','size':16}#, 'serif': ['computer modern roman']}
plt.rc('font',**font)

mpl.rcParams['lines.linestyle'] = '-'
mpl.rcParams['lines.linewidth'] = 1
##############################################
##############################################
# Path to folder containing .npz files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

# pathresp = "aia_resp_110809.sav"
# path = "D://emcubes_110809//emcubes_np//"
# # path = "D://emcubes_no_335_110809//emcubes_no_335_np//"
# # dir_emcubes = savedir + "emcubes_no335//"
# # xregion = [350,651]
# # yregion = [305,746]
# xregion = [410,601]
# yregion = [410,611]
# savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2011//"

##############################################

pathresp = "aia_resp_140910.sav"
path = "D://emcubes_140910//emcubes_np//"
# xregion = [250,750]
# yregion = [350,600]
savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2014//"
xregion = [350,650]
yregion = [380,530]
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

filelist = os.listdir(path)

xstart  = xregion[0]
xend    = xregion[1]
ystart  = yregion[0]
yend    = yregion[1]

color_list = [
    "#000000",  # black, 0 (converges)
    "#984ea3",  # purple, 1 (unbounded)
    "#ffff00",  # yellow, 2 (no soln) 
    "#FF63C2",  # pink, 3 (func did not converge)
    "#576BFF",  # blue,  10 (negative coeff from simplex) 
    "#88DC85",  # green, 11 (nocount; too little data)+
    ]
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']

color_map = colors.ListedColormap(color_list)
sat_color = colors.ListedColormap("#e41a1c")

allwaves  = [94,131,171,193,211,335]

##############################################
##############################################
# Functions
##############################################
##############################################

def imgs_emcubes(contour_sat = False, color_sat = False, status = False):
    """Produces frames for a movie of the event, where the data plotted is the sum of the emission cubes over logT.
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

        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        fig, ax = plt.subplots(figsize=(emcube.shape[1]*2*px, emcube.shape[0]*2*px),
                frameon = False)
        ax.set_axis_off()

        ax.imshow(emcube, cmap = 'gray', interpolation = 'none', origin = 'lower',)

        if contour_sat == True:

            sat = np.copy((vars["satmap"])[:,ystart:yend,xstart:xend])
            
            contour94 = ax.contour(sat[0,:,:], antialiased = False, colors = color_list[0])
            contour131 = ax.contour(sat[1,:,:], antialiased = False, colors = color_list[1])
            contour171 = ax.contour(sat[2,:,:], antialiased = False, colors = color_list[2])
            contour193 = ax.contour(sat[3,:,:], antialiased = False, colors = color_list[3])
            contour211 = ax.contour(sat[4,:,:], antialiased = False, colors = color_list[4])
            contour335 = ax.contour(sat[5,:,:], antialiased = False, colors = color_list[5])

            h0,_ = contour94.legend_elements()
            h1,_ = contour131.legend_elements()
            h2,_ = contour171.legend_elements()
            h3,_ = contour193.legend_elements()
            h4,_ = contour211.legend_elements()
            h5,_ = contour335.legend_elements()

            leg = ax.legend([h0[0], h1[0], h2[0], h3[0], h4[0], h5[0]],
                       ['$94 \\mathrm{\AA}$',
                    '$131 \\mathrm{\AA}$',
                    '$171 \\mathrm{\AA}$',
                    '$193 \\mathrm{\AA}$',
                    '$211 \\mathrm{\AA}$',
                    '$335 \\mathrm{\AA}$',],
                    # fontsize = 'x-large',
                    framealpha = 1)     
        
            for line in leg.get_lines():
                line.set_linestyle('-')

        if color_sat == True:
            satmap = np.sum((vars["satmap"])[:,ystart:yend, xstart:xend], axis = 0)

            satmasked = np.ma.masked_where(satmap < 1, satmap)      # makes array displaying which pixels are saturated
            ax.imshow(satmasked, cmap = sat_color, origin = 'lower', interpolation = 'none')

        if status == True: 
            # statuscube = vars.statuscube[ystart:yend, xstart:xend]
            statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]

            statusmasked = np.ma.masked_where(emcube != 0, statuscube)
            ax.imshow(statusmasked, cmap = color_map, vmin = 0, vmax = 5, origin = 'lower', interpolation = 'none')

        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        fig.savefig(dir_emcubes + file[:-4] + " " + str(i) + ".png", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    )        
        
        plt.clf()

    # make a test image to ensure the colors are as expected
    fig = plt.figure(num = 1, clear = True, frameon = False)
    fig.set_size_inches(emcube.shape[1], emcube.shape[0])
    ax = plt.Axes(fig, [0., 0., 1., 1.])
    ax.set_axis_off()
    fig.add_axes(ax)

    test_satarr = np.full_like(emcube, 0)
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
    
    plt.close()
    return

def imgs_statuscubes(contour_sat = False):
    """Produces frames for a movie of the event, where the data plotted a color-coded map of the error code
    returned by the DEM solver.
    Args:
        None.
    Returns:
        None. Saves a series of images into the directory savedir/statuscubes.
    """

    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file
        print(filename)

        vars = np.load(filename)

        statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]

        # 1: The objective function is unbounded.
        # 2: No solution satisfies the given constraints.
        # 3: The routine did not converge.
        statuscube[statuscube == 10] = 4    # SIMPLEX very occasionally returns negative coeffs even though positivity is a constraint. Now these are set to zero and the corresponding status code is 10.
        statuscube[statuscube == 11] = 5    # NOCOUNTS = WHERE(total(image,3) LT 10*eps). Set the status value of nocount pixels to 11.0

        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        fig, ax = plt.subplots(figsize=(statuscube.shape[1]*px, statuscube.shape[0]*px), frameon = False)
        ax.set_axis_off()
        
        ax.imshow(statuscube, cmap = color_map, vmin = 0, vmax = 5, interpolation = 'none', origin = 'lower')

        if contour_sat == True:
            sat = np.sum(np.copy((vars["satmap"])[:,ystart:yend,xstart:xend]), 0)
            sat = sat.astype('int32')
            ax.contour(sat, levels = [0, 2, 4], colors = contourclrs, antialiased = False)#linewidths = 100)

        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        fig.savefig(dir_statuscubes + file[:-4] + " " + str(i) + ".png", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    )        

        plt.clf()

    # make test image to ensure color maps work

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


def imgs_allchannels():
    """Produces frames for a movie of the event in all wavelengths.
    Args:
        None.
    Returns:
        None.
        Saves a series of images into the directory savedir/data_[wave]/.
    """

    colors = [] # list of colormaps
    for e in range(6):
        colors.append(aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom))

    # Calculate the correct figure size based on data dimensions
    data_width = xend - xstart
    data_height = yend - ystart
    height = [data_height, data_height]
    width = [data_width, data_width, data_width]

    print(data_width, data_height)
    print(data_width*3, data_height*2)

    px = 1/plt.rcParams['figure.dpi']  # pixel in inches

    fig, axs = plt.subplots(2, 3,
            figsize=(sum(width)*px, sum(height)*px),
            gridspec_kw={'height_ratios':height})

    titles = ['$94 \mathrm{\AA}$', 
            '$131 \mathrm{\AA}$', 
            '$171 \mathrm{\AA}$', 
            '$193 \mathrm{\AA}$', 
            '$211 \mathrm{\AA}$',
            '$335 \mathrm{\AA}$']

    axs = np.reshape(axs, (6))

    for i in range(len(filelist)):
        file = filelist[i]
        data = np.load(path + file)
        data = np.clip((data["datacube"])[:, ystart:yend, xstart:xend], a_min=0, a_max=None)**0.25

        print(i)
        for j in range(6):
            axs[j].imshow(data[j],
                        interpolation = 'none',
                        origin = 'lower',
                        cmap = colors[j])
            axs[j].axis('off')
            axs[j].text(0.02, 0.9,  # x, y
                        titles[j],  # string
                        transform = axs[j].transAxes,   # places on each axis individually
                        c = 'white',
                        size = 'xx-large',
                        )

        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        fig.savefig(savedir + "data_allchannels//" + file[:-4] + ".png")


def imgs_rel_err():

    for i in range(len(filelist)):

        file = filelist[i]
        filename = path + file
        print(filename)

        vars = np.load(filename)
        emcube = (np.sum( (vars["emcube"])[:,ystart:yend,xstart:xend], axis = 0)**0.25)

        px = 1/plt.rcParams['figure.dpi']  # pixel in inches
        fig, ax = plt.subplots(figsize=(emcube.shape[1]*2*px, emcube.shape[0]*2*px),
                frameon = False)
        ax.set_axis_off()

        ax.imshow(emcube, cmap = 'gray', interpolation = 'none', origin = 'lower',)

        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1, bottom=0, top=1)
        fig.savefig(dir_emcubes + file[:-4] + " " + str(i) + ".png", 
                    pad_inches = 0, 
                    bbox_inches= 'tight',
                    )        
        
        plt.clf()

    return