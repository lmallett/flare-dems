# (c) Lucien Mallett 2024

import scipy as sp
import numpy as np
import os

import matplotlib
# matplotlib.use("Agg")
matplotlib.use('TkAgg')

from matplotlib import rc
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

import matplotlib.pyplot as plt

from astropy.time import Time
from astropy.timeseries import TimeSeries
from astropy import visualization

import idl_colorbars
import astropy.units as u
import sunpy.visualization.colormaps.color_tables as aiacolormaps

import sunkit_image.trace as trace

import sunpy.io.special.genx

#####################################################
# Path to folder containing .sav files which contain image data

# OLD LOCATIONS
# path = "E:\\2011_08_09\\emcubes\\all\\"
# path = "E:\\2014_09_10\\emcubes\\all\\"
# path = "E:\\emcubes_110809\\uncropped\\"
# path = "E:\\emcubes_140910\\uncropped\\"

# path = "E:\\emcubes_110809\\"
# xregion = [350,650]
# yregion = [275,775]
# savedir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//"


path = "E:\\emcubes_140910\\"
xregion = [250,750]
yregion = [350,750]   
savedir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2014//"

# Path to response functions 
pathresp = "aia_resp_full.sav"

#####################################################

xstart  = xregion[0]
xend    = xregion[1]-1
ystart  = yregion[0]
yend    = yregion[1]-1

respinfo = sp.io.readsav(pathresp, python_dict = True, verbose = False)

logte = (((respinfo['r'])['A94'])[0])['logte'][0]
r94 = (((respinfo['r'])['A94'])[0])['tresp'][0]
r131 = (((respinfo['r'])['A131'])[0])['tresp'][0]
r171 = (((respinfo['r'])['A171'])[0])['tresp'][0]
r193 = (((respinfo['r'])['A193'])[0])['tresp'][0]
r211 = (((respinfo['r'])['A211'])[0])['tresp'][0]
r335 = (((respinfo['r'])['A335'])[0])['tresp'][0]

greenwhite = idl_colorbars.getcmap(8)
contourclrs = ['#ff14ec', '#8f0ad1', '#1e00b6']
filelist = os.listdir(path)

allwaves  = [94,131,171,193,211,335]

######################################################################
######################################################################
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

        datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
        # satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
        # satmasked = np.ma.masked_where(satmap < 1, satmap) # makes transparent satmap for overplotting
        emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])

        sel_ind = np.searchsorted(logte, vars.logt)
        logt = logte[sel_ind]

        print(emcube.shape)
        print(datacube.shape)

        # show image with line on pixel
        axs[1].imshow(np.sum(emcube, axis = 0)**0.25, origin = 'lower', cmap = 'binary_r')
        # axs[1].scatter(x,y, c = 'aqua', marker = '+')

        axs[1].axhline(y, color="r", linestyle="--")
        axs[1].axvline(x, color="r", linestyle="--")
        axs[1].axis('off')

        # IDL is col-major, so the data is unfortunately stored as [temp, y, x]
        axs[0].step(logt, emcube[:,y,x], where = 'post', color = 'k', label = 'EM')

        axs[0].tick_params(axis='both', which = 'both',
                           direction = 'in',
                           bottom = True, top= True, left= True, right = True)
        axs[0].tick_params(which = 'major', length = 10)
        axs[0].tick_params(which = 'minor', length = 5)

        axs[0].set_xlim(5.6, 8.0)
        axs[0].set_ylim(1, 10**6)
        axs[0].set_xlabel("$log(T)$")
        axs[0].set_ylabel("Emission Measure [$10^{26} \mathrm{cm}^{-5}$]")
        axs[0].set_yscale("log")
        axs[0].set_aspect(0.2)
        axs[0].set_xticks(np.arange(5.6, 8.0, 0.2))

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
        plt.savefig(fname = savedir + "lociplots_280-140//" + file[:-4] + ".png")#, transparent = True)

        # fig.show()


######################################################################
######################################################################
######################################################################
######################################################################

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
    fig.set_size_inches(11, 5.5)
    fig.tight_layout()
    plt.subplots_adjust(left = 0.07)

    for i in range(len(filelist)):

        axs.clear()

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        print(file)

        vars = sp.io.readsav(filename)

        datacube = np.copy(vars.datacube[:,ystart:yend, xstart:xend])
        satmap = np.copy(np.sum(vars.satmap[:,ystart:yend, xstart:xend], axis=0))
        emcube = np.copy(vars.emcube[:,ystart:yend, xstart:xend])

        sel_ind = np.searchsorted(logte, vars.logt)
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
        axs.set(title = f"num. saturated channels: {int(satmap[y,x])}")

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
        axs.legend(fontsize = 20, loc = 'upper left')
        plt.pause(0.001)        
        # plt.savefig(fname = savedir + "lociplots_only_144-214//" + file[:-4] + ".png")#, transparent = True)
        # plt.savefig(fname = savedir + "lociplots_only_153-214//" + file[:-4] + ".png")#, transparent = True)
        plt.savefig(fname = savedir + "lociplots_only_280_140//" + file[:-4] + ".png")#, transparent = True)


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
    # plt.style.use("dark_background")
    plt.savefig(fname = "C:\\Users\\Lucien\\Documents\\School\\Research - Thesis\\response.png")#, transparent = True)

######################################################################
######################################################################
######################################################################
######################################################################
    
def exptimes(wavelength = allwaves):
    """Produces a graph of the exposure times for each channel over time.

    Keyword Arguments:
        wavelength -- An integer or list of integers of channel wavelengths in Angstroms. (default: {allwaves}, all available channels.)
    """

    if type(wavelength) != list:
        wavelength = [wavelength]

    n = len(filelist)

    exptimes = np.zeros([n, 6])
    times = []
    print("done")

    for i in range(n):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        time = file[11:-4].replace("_",":")
        times.append(time)

        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)
        exptimes[i] = np.copy(vars.exptimeout)
        print(filename)

    ts = TimeSeries(time = times)
    ts.sort('time')
    ts['times_94'] = exptimes[:,0]
    ts['times_131'] = exptimes[:,1]
    ts['times_171'] = exptimes[:,2]
    ts['times_193'] = exptimes[:,3]
    ts['times_211'] = exptimes[:,4]
    ts['times_335'] = exptimes[:,5]

    print(ts)

    print(ts.time.format)
    fig, ax = plt.subplots()

    with visualization.time_support(format = 'isot', simplify=True):  

        ax.plot(ts.time, ts['times_94'], label = '94')
        ax.plot(ts.time, ts['times_131'], label = '131')
        ax.plot(ts.time, ts['times_171'], label = '171')
        ax.plot(ts.time, ts['times_193'], label = '193')
        ax.plot(ts.time, ts['times_211'], label = '211')
        ax.plot(ts.time, ts['times_335'], label = '335')

        # ax.xlabel('Julian Date')
        # ax.ylabel('SAP Flux (e-/s)')
        ax.grid()
    plt.legend()
    
    plt.show()

######################################################################
######################################################################
######################################################################
######################################################################

def num_sat():
    """Produces a "number saturated vs. time" graph for a given data set.
    """

    n = len(filelist)

    numSat = np.zeros([n, 6])
    norms = np.zeros([n, 6])

    times = []

    for i in range(n):

        # read in data & modify it
        file = filelist[i]
        filename = path + file

        time = file[11:-4].replace("_",":")
        times.append(time)

        vars = sp.io.readsav(filename)#, python_dict = False, verbose = False)
        sums = np.sum(vars.satmap, axis = (1,2))

        numSat[i] = sums

        # Catch for when there are no saturated channels by setting them to 0
        norms[i] = np.nan_to_num(sums / np.linalg.norm(sums))

        print(sums)

    ts = TimeSeries(time = times)
    ts['times_94'] = numSat[:,0]
    ts['times_131'] = numSat[:,1]
    ts['times_171'] = numSat[:,2]
    ts['times_193'] = numSat[:,3]
    ts['times_211'] = numSat[:,4]
    ts['times_335'] = numSat[:,5]

    ts['norm_94'] = norms[:,0]
    ts['norm_131'] = norms[:,1]
    ts['norm_171'] = norms[:,2]
    ts['norm_193'] = norms[:,3]
    ts['norm_211'] = norms[:,4]
    ts['norm_335'] = norms[:,5]

    ts.sort('time')
    print(ts)
    
    # print("Average times when there are saturated pixels:")
    # print("94A:" np.average())

    fig, ax = plt.subplots(2, sharex = True)

    with visualization.time_support(format = 'isot', simplify=True):  

        ax[0].plot(ts.time, ts['times_94'], label = "$94 \\mathrm{\AA}$", linestyle = 'dashed', c = 'r')
        ax[0].plot(ts.time, ts['times_131'], label = "$131 \\mathrm{\AA}$", c = 'y')
        ax[0].plot(ts.time, ts['times_171'], label = "$171 \\mathrm{\AA}$", linestyle = 'dashed', c = 'g')
        ax[0].plot(ts.time, ts['times_193'], label = "$193 \\mathrm{\AA}$", c = 'c')
        ax[0].plot(ts.time, ts['times_211'], label = "$211 \\mathrm{\AA}$", linestyle = 'dashed', c = 'm')
        ax[0].plot(ts.time, ts['times_335'], label = "$335 \\mathrm{\AA}$", c = 'k')
        ax[0].grid()
        ax[0].legend()

        ax[1].plot(ts.time, ts['norm_94'], label = "$94 \\mathrm{\AA}$", linestyle = 'dashed', c = 'r')
        ax[1].plot(ts.time, ts['norm_131'], label = "$131 \\mathrm{\AA}$", c = 'y')
        ax[1].plot(ts.time, ts['norm_171'], label = "$171 \\mathrm{\AA}$", linestyle = 'dashed', c = 'g')
        ax[1].plot(ts.time, ts['norm_193'], label = "$193 \\mathrm{\AA}$", c = 'c')
        ax[1].plot(ts.time, ts['norm_211'], label = "$211 \\mathrm{\AA}$", linestyle = 'dashed', c = 'm')
        ax[1].plot(ts.time, ts['norm_335'], label = "$335 \\mathrm{\AA}$", c = 'k')
        ax[1].grid()
        ax[1].legend()
        ax[1].set(ylim = [-0.1,1.1])

    plt.show()


######################################################################
######################################################################
######################################################################
######################################################################

# We tried to use this, initially, to detect loop locations. However, this was not consistent enough for our purposes.

# directories to put individual frames in
dirs = {
    94: "data_94/",
    131: "data_131/",
    171: "data_171/",
    193: "data_193/",
    211: "data_211/",
    335: "data_335/",
    "emcube": "emcubes//"}

def looptrace(wavelength = allwaves):
    """Produces frames for a movie of the event in the specified AIA wavelength (94, 131, 171, 193, 211, or 335), with lines
       shown in the locations where loops are detected.
    Args:
        (Optional) An integer, or list of integers, that is a subset of [94,131,171,193,211,335].
    Returns:
        None.
        Saves a series of images into the directory savedir/data_[wave]/.
    """

    if type(wavelength) != list:
        wavelength = [wavelength]

    imglist = []    # list of datacubes
    for i in range(25):#len(filelist)):

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