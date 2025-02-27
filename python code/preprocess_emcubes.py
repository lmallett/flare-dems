import os
import scipy as sp
import numpy as np

# path = "D://emcubes_110809//emcubes//"
# savdir = "C://Users//Lucien//Documents//School//Research - Thesis//emcubes_np//"

# path = "D://emcubes_no_335_110809//emcubes_no_335//"
# savdir = "D://emcubes_no_335_110809//emcubes_no_335_np//"

# path = "D://emcubes_110809//emcubes_mchannel//m335//"
# savdir = "D://emcubes_110809//emcubes_mchannel//m335_np//"

path = "D://emcubes_140910//emcubes//"
savdir = "D://emcubes_140910//emcubes_np//"
# path = "D://emcubes_140910//emcubes_no_335//"
# savdir = "D://emcubes_140910//emcubes_no_335_np//"

filelist = os.listdir(path)

for i in range(len(filelist)):

    file = filelist[i]
    vars = sp.io.readsav(path + file)
    print(file)

    emcube      = vars.emcube
    statuscube  = vars.statuscube
    logt        = vars.logt
    datacube    = vars.datacube
    waveout     = vars.waveout
    exptimeout  = vars.exptimeout
    satmap      = vars.satmap

    np.savez(savdir + file[:-4],
             emcube = emcube,
             statuscube = statuscube,
             logt = logt,
             datacube = datacube,
             waveout = waveout,
             exptimeout = exptimeout,
             satmap = satmap
             )