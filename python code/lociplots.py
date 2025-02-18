# (c) Lucien Mallett 2024

import scipy as sp
import numpy as np
import os

import matplotlib
# matplotlib.use("Agg")
matplotlib.use('TkAgg')

# from matplotlib import rc
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import matplotlib.pyplot as plt
import idl_colorbars

from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy import visualization
import astropy.units as u

import sunpy.visualization.colormaps.color_tables as aiacolormaps
import sunkit_image.trace as trace
import sunpy.io.special.genx

##############################################
##############################################
# Path to folder containing .sav files, which contain image data at different timesteps
# Select region for view
##############################################
##############################################

pathresp = "aia_resp_110809.sav"
path = "D://emcubes_110809//emcubes//"
# path = "D://emcubes_no_335_110809//emcubes_no_335//"
xregion = [350,650]
yregion = [275,775]
savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2011//"

pathresp = "aia_resp_140910.sav"
path = "D://emcubes_140910//emcubes_np//"
# path = "D://emcubes_140910//emcubes_no_335//"
# xregion = [250,750]
# yregion = [350,600]
xregion = [375,575]
yregion = [400,480]
savedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2014//"

######################################################################
######################################################################
# Initialization
######################################################################
######################################################################

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

# load in response functions as 20x1 arrays
respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

# logte = (((respinfo['r'])['A94'])[0])['logte'][0]
# r94 = (((respinfo['r'])['A94'])[0])['tresp'][0]
# r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
# r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
# r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
# r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
# r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]
logte= (respinfo["aia_tresp"][0])[0]
r94  = (respinfo["aia_tresp"][0])[2]
r131 = (respinfo["aia_tresp"][0])[3]
r171 = (respinfo["aia_tresp"][0])[4]
r193 = (respinfo["aia_tresp"][0])[5]
r211 = (respinfo["aia_tresp"][0])[6]
r335 = (respinfo["aia_tresp"][0])[7]

allwaves  = [94,131,171,193,211,335]

greenwhite = idl_colorbars.getcmap(8)
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
filelist = os.listdir(path)

######################################################################
######################################################################
# Plotting
######################################################################
######################################################################  

def lociplots(x, y):
    """Produces images saved in savedir/(specify by editing code) of the lociplots at pixel x, y in the image, one for each timestep.
    The images also show the location of the pixel by displaying the solar image on the right.

    Arguments:
        x -- x pixel specified.
        y -- y pixel specified.
    """

    global sel_ind

    resp = [r94, r131, r171, r193, r211, r335]

    linecolors = ['r','g','b', 'c', 'm', 'y']
    fig, axs = plt.subplots(1,2)
    fig.set_size_inches(11, 4)
    fig.tight_layout()
    plt.subplots_adjust(left = 0.07)

    for i in range(len(filelist)):

        axs[0].clear()

        # read in data & modify it
        file = filelist[i]

        print(file)

        if file[-4:] == ".sav":
            vars = sp.io.readsav(path + file)
            datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
            # satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
            # satmasked = np.ma.masked_where(satmap < 1, satmap) # makes transparent satmap for overplotting
            emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])
            sel_ind = np.searchsorted(logte, vars.logt)

        if file[-4:] == ".npz":
            vars = np.load(path + file)
            datacube = (vars["datacube"])[:,ystart:yend, xstart:xend]
            emcube = (vars["emcube"])[:,ystart:yend, xstart:xend]
            # satmap = np.sum((vars["satmap"])[:, ystart:yend, xstart:xend], axis=0)
            # statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]
            # statuscube = np.ma.masked_where(
            #     satmap != 0, statuscube
            # )  # exclude pixels that were saturated
            sel_ind = np.searchsorted(logte, vars["logt"])

        logt = logte[sel_ind]

        print(emcube.shape)
        print(datacube.shape)

        # show image with line on pixel
        axs[1].imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = 'binary_r')
        # axs[1].scatter(x,y, c = 'aqua', marker = '+')

        # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
        axs[0].step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

        axs[0].set(aspect = 0.2,
                title = "Loci curves, data $y_i$ divided by response function $K_i$",
                xlabel = r"$\mathrm{log(T)}$",
                ylabel = "Emission Measure $[10^{26} \;\mathrm{cm}^{-5}]$",
                xlim = [5.6, 8.0],
                ylim = [1, 1e6],
                xticks = np.arange(5.6, 8.0, 0.2),
                yscale = "log")
        axs[0].tick_params(axis='both', which = 'both', direction = 'in',
                           bottom = True, top= True, left= True, right = True)
        axs[0].tick_params(which = 'major', length = 10)
        axs[0].tick_params(which = 'minor', length = 5)
        axs[1].axhline(y, color="r", linestyle="--")
        axs[1].axvline(x, color="r", linestyle="--")
        axs[1].axis('off')

        for j in range(6):

            respj = resp[j]
            respj = respj[sel_ind] # 20x1 array

            # plot responses
            # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
            datarray = np.full_like(respj, datacube[j,y,x])
            datarray2 = datarray / respj
            datarray2 = ( datarray2/(10**(26)) ) # 20x1 array

            axs[0].semilogy(logt, datarray2, color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))

        # plot emission cube
        axs[0].legend(loc = 'upper left')
        plt.pause(0.001)
        
        # plt.style.use("dark_background")
        # plt.savefig(fname = savedir + "lociplots_144-214//" + file[:-4] + ".png")#, transparent = True)
        # plt.savefig(fname = savedir + "lociplots_153-214//" + file[:-4] + ".png")#, transparent = True)
        # plt.savefig(fname = savedir + "lociplots_260-75//" + file[:-4] + ".png")#, transparent = True)
        plt.savefig(fname = savedir + "lociplots_125-25//" + file[:-4] + ".png", dpi = 300)#, transparent = True)

        # fig.show()


######################################################################
######################################################################
######################################################################
######################################################################


def imgs_only(x, y):
    """Produces images saved in savedir/(specify by editing code) of the solar image at pixel x, y in the image, one for each timestep.
    Arguments:
        x -- x pixel specified.
        y -- y pixel specified.
    """

    global sel_ind
    fig, axs = plt.subplots()
    fig.set_size_inches(6, 3)
    plt.tight_layout()
    plt.subplots_adjust(left = 0.05,
                        right = 0.95,
                        bottom = 0.05,
                        top = 0.95)

    for i in range(len(filelist)):

        axs.clear()

        # read in data & modify it
        file = filelist[i]

        print(file)

        if file[-4:] == ".sav":
            vars = sp.io.readsav(path + file)
            emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])
            sel_ind = np.searchsorted(logte, vars.logt)

        if file[-4:] == ".npz":
            vars = np.load(path + file)
            emcube = (vars["emcube"])[:,ystart:yend, xstart:xend]
            sel_ind = np.searchsorted(logte, vars["logt"])

        # show image with line on pixel
        axs.imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = 'binary_r')
        axs.axhline(y, color="r", linestyle="--")
        axs.axvline(x, color="r", linestyle="--")
        axs.axis('off')
        
        plt.savefig(fname = savedir + "lociplots_125-25_imgonly//" + file[:-4] + ".png", dpi = 300)#, bbox_inches="tight")

def lociplots_only(x, y):
    """Produces images saved in savedir/(specify by editing code) of the lociplots at pixel x, y in the image, one for each timestep.

    Arguments:
        x -- x pixel specified.
        y -- y pixel specified.
    """

    global sel_ind

    resp = [r94, r131, r171, r193, r211, r335]

    linecolors = ['r','g','b', 'c', 'm', 'y']
    fig, axs = plt.subplots()
    fig.set_size_inches(6, 3)
    fig.tight_layout()
    plt.subplots_adjust(left = 0.1)

    for i in range(len(filelist)):

        axs.clear()

        # read in data & modify it
        file = filelist[i]
        print(file)
        if file[-4:] == ".sav":
            vars = sp.io.readsav(path + file)
            datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
            emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])
            # satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
            sel_ind = np.searchsorted(logte, vars.logt)

        if file[-4:] == ".npz":
            vars = np.load(path + file)
            datacube = (vars["datacube"])[:,ystart:yend, xstart:xend]
            emcube = (vars["emcube"])[:,ystart:yend, xstart:xend]
            satmap = np.sum((vars["satmap"])[:, ystart:yend, xstart:xend], axis=0)
            # statuscube = (vars["statuscube"])[ystart:yend, xstart:xend]
            # statuscube = np.ma.masked_where(
            #     satmap != 0, statuscube
            # )  # exclude pixels that were saturated
            sel_ind = np.searchsorted(logte, vars["logt"])

        logt = logte[sel_ind]

        # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
        axs.step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

        axs.tick_params(axis='both', which = 'both',
                           direction = 'in',
                           bottom = True, top= True, left= True, right = True)
        axs.tick_params(which = 'major', length = 10)
        axs.tick_params(which = 'minor', length = 5)

        axs.set_xlim(5.6, 8.0)
        axs.set_ylim(1, 10**6)
        axs.set_xlabel("$log(T)$")
        axs.set_ylabel("Emission Measure [$10^{26} \mathrm{cm}^{-5}$]")
        axs.set_yscale("log")
        axs.set_aspect(0.2)
        axs.set_xticks(np.arange(5.6, 8.0, 0.2))
        # axs.set(title = f"num. saturated channels: {int(satmap[y,x])}")

        for j in range(6):

            respj = resp[j]
            respj = respj[sel_ind] # 20x1 array

            # plot responses
            # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
            datarray = np.full_like(respj, datacube[j,y,x])
            datarray2 = datarray / respj
            datarray2 = ( datarray2/(10**(26)) ) # 20x1 array

            axs.semilogy(logt, datarray2, color = linecolors[j], label = "${0} {1}$".format(allwaves[j], "\mathrm{\AA}"))

        # plot
        axs.legend(loc = 'upper left')
        plt.pause(0.001)        
        # plt.savefig(fname = savedir + "lociplots_only_144-214//" + file[:-4] + ".png")#, transparent = True)
        # plt.savefig(fname = savedir + "lociplots_only_153-214//" + file[:-4] + ".png")#, transparent = True)
        plt.savefig(fname = savedir + "lociplots_only_125-25//" + file[:-4] + ".png", dpi = 300)#, transparent = True)


######################################################################
######################################################################
######################################################################
######################################################################

# We tried to use this, initially, to detect loop locations. However, this was not consistent enough for our purposes.

# directories to put individual frames in
# dirs = {
#     94: "data_94/",
#     131: "data_131/",
#     171: "data_171/",
#     193: "data_193/",
#     211: "data_211/",
#     335: "data_335/",
#     "emcube": "emcubes//"}

# def looptrace(wavelength = allwaves):
#     """Produces frames for a movie of the event in the specified AIA wavelength (94, 131, 171, 193, 211, or 335), with lines
#        shown in the locations where loops are detected.
#     Args:
#         (Optional) An integer, or list of integers, that is a subset of [94,131,171,193,211,335].
#     Returns:
#         None.
#         Saves a series of images into the directory savedir/data_[wave]/.
#     """

#     if type(wavelength) != list:
#         wavelength = [wavelength]

#     imglist = []    # list of datacubes
#     for i in range(25):#len(filelist)):

#         # read in data & modify it
#         file = filelist[i]
#         filename = path + file

#         print(file)
#         vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)

#         data = np.copy(vars.datacube)
#         data[data < 0] = 0
#         # data = data**0.25

#         # save modifed datacubes
#         imglist.append(data)

#     # maps what wavelength we want to 0-5
#     # convenient for accessing datacubes
#     indices = [2]#allwaves.index(wave) for wave in wavelength]

#     for e in indices:

#         colors = aiacolormaps.aia_color_table(allwaves[e] * u.Angstrom) 
#         print(e)

#         for i in range(len(filelist)):            

#             # get data at timestamp, for only the relevant channel
#             file = filelist[i]
#             data = (imglist[i])[e,:,:]

#             # loops = trace.occult2(data, nsm1=3, rmin=30, lmin=25, nstruc=1000, ngap=0, qthresh1=0.0, qthresh2=3.0)
#             loops = trace.occult2(data, nsm1=3, rmin=10, lmin=25, nstruc=1000, ngap=0, qthresh1=0.0, qthresh2=3.0)

#             print('loop done')
#             fig = plt.figure(frameon = False)
#             fig.set_size_inches(data.shape[1], data.shape[0])
#             ax = plt.Axes(fig, [0., 0., 1., 1.])
#             ax.set_axis_off()
#             fig.add_axes(ax)

#             data = data**0.25
#             ax.imshow(data, cmap = colors, interpolation='none', origin='lower')

#             print('about to plot loops')
#             for loop in loops:
#                 loop = np.array(loop)
#                 x, y = loop.T
#                 plt.scatter(x, y, c = 'blue', s = 10000)

#             print('saving')
#             # i have no idea why mpl is so picky about fig vs plt and so on. changing stuff made me use 22GB of memory before i killed it.
#             fig.savefig(savedir + "loops_" + dirs[allwaves[e]] + file[:-4] + ".png", dpi = 1)

#     plt.close()