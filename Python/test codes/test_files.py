import os
import numpy as np
import scipy as sp

path = "E:\\2011_08_09.tar\\2011_08_09\\emcubes\\all\\"

filelist = os.listdir(path)
file = filelist[15]
filename = path + file

vars = sp.io.readsav(filename)
data = np.copy(vars.datacube)
# wheredata = np.argwhere(data > 15000)

# print(data[wheredata])
print(vars.datacube)
# print(vars.satmap)
