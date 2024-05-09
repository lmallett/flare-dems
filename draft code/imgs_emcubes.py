import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import idl_colorbars
import os

# Path to folder containing .sav files, which contain image data
path = "E:\\2011_08_09\emcubes\\all\\"

greenwhite = idl_colorbars.getcmap(8)
filelist = os.listdir(path)

fig, ax = plt.subplots()

imglist = []

for i in range(len(filelist)):
    file = filelist[i]
    filename = path + file

    print(file)
    vars = sp.io.readsav(path + file, python_dict = False, verbose = False)
    # print(vars.keys())

    data = np.sum(vars.emcube, axis = 0)**0.25
    # np.info(a)

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
ani.save("movie.mp4")
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


