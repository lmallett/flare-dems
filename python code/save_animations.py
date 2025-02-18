# Lucien Mallett (c)
# Produces animated images from still frames.

import imageio
from PIL import Image
import os

import cv2
import os

basedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2014//"
subdirs = ["data_94//", "data_131//", "data_171//", "data_193//", "data_211//", "data_335//", "emcubes//", "statuscubes//"]

def save_gif(dir, name):
    """Produces a GIF of the images in a folder.
    Input: the directory containing the images to be saved."""
    filelist = [ f for f in os.listdir(dir) if f.endswith('.png') ]
    imglist = []
    for filename in filelist:
        print(filename)
        imglist.append(imageio.imread(dir + filename))
    imageio.mimsave(dir + str(name) + ".gif", imglist)

def save_webp(dir, name):
    """Produces an animated .webp of the images in a folder.
    Input: the directory containing the images to be saved."""

    filelist = [ f for f in os.listdir(dir) if f.endswith('.png') ]
    imglist = []
    
    for filename in filelist:
        im = Image.open(dir + filename)
        imglist.append(im)

    imglist[0].save(dir +  str(name) + ".webp",
        'webp',
        save_all=True, 
        lossless= True,
        quality = 100,
        append_images = imglist[1:], 
        optimize= False,
        duration = 100,
        loop    = 0)


def save_mp4(dir, name):
    """Produces an animated .mp4 of the images in a folder.
    Input: the directory containing the images to be saved."""
    filelist = [ f for f in os.listdir(dir) if f.endswith('.png') ]

    # Read the first image to get the dimensions
    frame = cv2.imread(os.path.join(dir, filelist[0]))
    height, width, _ = frame.shape

    # Create a VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Use a suitable codec
    video = cv2.VideoWriter(dir + str(name) + ".mp4", fourcc, 10, (width, height))  # 20 is the frame rate
    cv2.VIDEOWRITER_PROP_QUALITY = 100

    # Write each image to the video
    for image in filelist:
        frame = cv2.imread(os.path.join(dir, image))
        video.write(frame)

    print(f"Total images to process: {len(filelist)}")
    print(f"Output video name: {name}")

    # Release the video writer
    video.release()

for subdir in subdirs:
    dir = basedir + subdir
    save_mp4(dir, subdir[:-2])

# Alternate attempt

# import os
# import subprocess

# # Define directory and video name
# basedir = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//2014//statuscubes//"
# name = "output_video"

# filelist = [f for f in os.listdir(basedir) if f.endswith('.png')]

# # Create a temporary file list for ffmpeg
# filelist_path = os.path.join(basedir, "filelist.txt")
# with open(filelist_path, "w") as file:
#     for img in filelist:
#         file.write(f"file '{os.path.join(basedir, img)}'\n")

# # Run ffmpeg to create a lossless H.264 video
# output_video = os.path.join(basedir, f"{name}.mp4")
# ffmpeg_command = [
#     'ffmpeg',
#     '-f', 'concat', '-safe', '0', '-i', filelist_path,
#     '-c:v', 'libx264', '-preset', 'ultrafast', '-crf', '0', '-pix_fmt', 'yuv420p',  # Lossless mode with H.264
#     output_video
# ]

# # Execute the ffmpeg command
# subprocess.run(ffmpeg_command)

# # Remove the temporary file list
# os.remove(filelist_path)

# print(f"Total images to process: {len(filelist)}")
# print(f"Output video name: {output_video}")