import scipy as sp
import numpy as np
import os

import matplotlib
# matplotlib.use("Agg")
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

import idl_colorbars
import astropy.units as u
import sunpy.visualization.colormaps.color_tables as aiacolormaps

import sunkit_image.trace as trace

import sunpy.io.special.genx

# Path to folder containing .sav files, which contain image data

path = "E:\\2011_08_09\\emcubes\\all\\"
path = "E:\\2014_09_10\\emcubes\\all\\"


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


# def lociplot(coords = [125,175]):

def looptrace(wavelength = allwaves):
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
        # data = data**0.25

        # save modifed datacubes
        imglist.append(data)

    # maps what wavelength we want to 0-5
    # convenient for accessing datacubes
    indices = [2]#allwaves.index(wave) for wave in wavelength]

    for e in indices:

        colors = aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom) 
        print(e)

        for i in range(len(filelist)):            

            # get data at timestamp, for only the relevant channel
            file = filelist[i]
            data = (imglist[i])[e,:,:]

            # loops = trace.occult2(data, nsm1=3, rmin=30, lmin=25, nstruc=1000, ngap=0, qthresh1=0.0, qthresh2=3.0)
            loops = trace.occult2(data, nsm1=3, rmin=10, lmin=25, nstruc=1000, ngap=0, qthresh1=0.0, qthresh2=3.0)

            print('loop done')
            fig = plt.figure(frameon = False)
            fig.set_size_inches(data.shape[1], data.shape[0])
            ax = plt.Axes(fig, [0., 0., 1., 1.])
            ax.set_axis_off()
            fig.add_axes(ax)

            data = data**0.25
            ax.imshow(data, cmap = colors, interpolation='none', origin='lower')

            print('about to plot loops')
            for loop in loops:
                loop = np.array(loop)
                x, y = loop.T
                plt.scatter(x, y, c = 'blue', s = 10000)

            print('saving')
            # i have no idea why mpl is so picky about fig vs plt and so on. changing stuff made me use 22GB of memory before i killed it.
            fig.savefig(savedir + "loops_" + dirs[allwaves[e]] + file[:-4] + ".png", dpi = 1)

    plt.close()

######################################################################
######################################################################
######################################################################
######################################################################

# select coronal loop pixels
# pick the brightest pixels in the loop that converge

# find loci plots of those
# pick non coronal loop pixels in the ROI
# find loci plots of those    

    # logt=indgen(20)*0.1+5.6
    # r=aia_get_response(/dn,/temperature,/full,/evenorm,timedepend_date='30-Jul-2021 18:00:00')
    # r94= r.a94.tresp
    # r131=r.a131.tresp
    # r171=r.a171.tresp
    # r193=r.a193.tresp
    # r211=r.a211.tresp
    # r335=r.a335.tresp

    # nch=6

    # sel_ind=value_locate(r.a94.logte,logt)        

def lociplots(x, y):

    # load in response functions as 20x1 arrays
    respinfo = sp.io.readsav("aia_response_2011.sav", python_dict = False, verbose = False)
    r94 = respinfo['r94']
    r131 = respinfo['r131']
    r171 = respinfo['r171']
    r193 = respinfo['r193']
    r211 = respinfo['r211']
    r335 = respinfo['r335']
    logt = respinfo['logt']
    logte = ((respinfo['r'])['logte'])[0]

    print(logte)

    # only select data at the same logT values (response functions are defined at more logT values than DEM solutions)
    sel_ind = respinfo['sel_ind']
    resp = [r94, r131, r171, r193, r211, r335]
    linecolors = ['k', 'r', 'g', 'b', 'c', 'pink', 'y']

    linecolors = ['r','g','b', 'c', 'm', 'y']
    fig, axs = plt.subplots(1,2)
    fig.set_size_inches(11, 5.5)
    fig.tight_layout()
    plt.subplots_adjust(left = 0.07)

    for i in range(len(filelist)):

        axs[0].clear()

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        print(file)

        vars = sp.io.readsav(filename)
        datacube = np.copy(vars.datacube)
        emcube = np.copy(vars.emcube)
        
        print(emcube.shape)
        print(datacube.shape)

        # show image with + on pixel
        axs[1].imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = greenwhite)
        axs[1].scatter(x,y, c = 'aqua', marker = '+')

        # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
        axs[0].step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

        axs[0].tick_params(axis='both', which = 'both',
                           direction = 'in',
                           bottom = True, top= True, left= True, right = True)
        axs[0].tick_params(which = 'major', length = 10)
        axs[0].tick_params(which = 'minor', length = 5)

        axs[0].set_xlim(5.6, 7.5)
        axs[0].set_ylim(1, 10**6)
        axs[0].set_xlabel("$log(T)$")
        axs[0].set_ylabel("Emission Measure [$10^{26} cm^{-5}$]")
        axs[0].set_yscale("log")
        axs[0].set_aspect(0.2)
        axs[0].set_xticks(np.arange(5.6, 7.5, 0.2))

        for j in range(6):

            respj = resp[j]
            respj = respj[sel_ind] # 20x1 array

            # plot responses
            # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
            datarray = np.full_like(respj, datacube[j,y,x])
            datarray2 = datarray / respj
            datarray2 = ( datarray2/(10**(26)) ) # 20x1 array

            axs[0].semilogy(logt, datarray2, color = linecolors[j], label = allwaves[j])

        # plot emission cube
        axs[0].legend(loc = 'upper left')
        plt.pause(0.001)
        plt.savefig(fname = savedir + "lociplots_140-250/" + file[:-4] + ".png")
        # fig.show()

######################################################################
######################################################################
######################################################################
######################################################################

def norm_temp_responses():
    """Plot the temperature response functions for AIA, normalized, to show where each channel is most sensitive."""

    # load in response functions as 101x1 arrays
    respinfo = sp.io.readsav("aia_response_2011.sav", python_dict = False, verbose = True)
    r94 = respinfo['r94']
    r131 = respinfo['r131']
    r171 = respinfo['r171']
    r193 = respinfo['r193']
    r211 = respinfo['r211']
    r335 = respinfo['r335']
    logte = ((respinfo['r'])['logte'])[0]

    #normalize the functions
    r94 = r94 / max(r94)
    r131 = r131 / max(r131)
    r171 = r171 / max(r171)
    r193 = r193 / max(r193)
    r211 = r211 / max(r211)
    r335 = r335 / max(r335)

    resp = [r94, r131, r171, r193, r211, r335]
    linecolors = ['r','g','b', 'c', 'm', 'y']

    fig, axs = plt.subplots()
    # # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
    # axs.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

    axs.minorticks_on()
    axs.tick_params(axis='both', which = 'both',
                        direction = 'in',
                        bottom = True, top= True, left= True, right = True)
    axs.tick_params(which = 'major', length = 10)
    axs.tick_params(which = 'minor', length = 5)

    axs.set_xlim(4.5, 7.5)
    axs.set_ylim(0, 1.1)
    axs.set_xlabel(r"$\log_{10}(T/\mathrm{K})$")
    axs.set_ylabel("Normalized AIA Temp. Response")

    for j in range(6):

        respj = resp[j]
        
        # plot responses
        axs.plot(logte, respj, color = linecolors[j], label = allwaves[j])
    
    axs.legend()
    plt.show()

######################################################################
######################################################################
######################################################################
######################################################################

def temp_responses():
    """Plot the temperature response functions for AIA to show where each channel is most sensitive."""

    # load in response functions as 101x1 arrays
    respinfo = sp.io.readsav("aia_response_2011.sav", python_dict = False, verbose = True)
    r94 = respinfo['r94']
    r131 = respinfo['r131']
    r171 = respinfo['r171']
    r193 = respinfo['r193']
    r211 = respinfo['r211']
    r335 = respinfo['r335']
    logte = ((respinfo['r'])['logte'])[0]

    resp = [r94, r131, r171, r193, r211, r335]
    linecolors = ['r','g','b', 'c', 'm', 'y']

    fig, axs = plt.subplots()

    axs.minorticks_on()
    axs.tick_params(axis='both', which = 'both',
                        direction = 'in',
                        bottom = True, top= True, left= True, right = True)
    axs.tick_params(which = 'major', length = 10)
    axs.tick_params(which = 'minor', length = 5)

    axs.set_xlim(4.5, 7.5)
    axs.set_ylim(10**(-28), 10**(-23))
    axs.set_xlabel(r"$\log_{10}(T/\mathrm{K})$")
    axs.set_ylabel(r"AIA Temp. Response $(\mathrm{DN}\;\mathrm{cm}^5 \; \mathrm{px}^{-1})$")

    for j in range(6):
        
        # plot responses
        respj = resp[j]
        axs.semilogy(logte, respj, color = linecolors[j], label = allwaves[j])
    
    axs.legend()
    plt.show()