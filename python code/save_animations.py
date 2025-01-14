# Lucien Mallett (c)
# Produces animated images from still frames.

import imageio
from PIL import Image
import os

dir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//emcubes//"
dir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//emcubes_no335//"
dir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//contoured_emcubes//"
# dir = "C://Users//Lucien//Documents//School//Research - Thesis//movies//2011//timing_tests//"

def save_gif():
    """Produces a GIF of the images in a folder.
    Input: the directory containing the images to be saved."""
    filelist = os.listdir(dir)
    imglist = []
    for filename in filelist:
        print(filename)
        imglist.append(imageio.imread(dir + filename))
    imageio.mimsave(dir + "emcubes.gif", imglist)

def save_webp():
    """Produces an animated .webp of the images in a folder.
    Input: the directory containing the images to be saved."""

    filelist = os.listdir(dir)
    imglist = []
    for file in filelist:
        im = Image.open(dir + file)
        imglist.append(im)
    imglist[0].save(dir + "emcube.webp",
        'webp',
        save_all=True, 
        lossless= True,
        quality = 100,
        append_images = imglist[1:], 
        optimize= False,
        duration = 100,
        loop    = 0)
