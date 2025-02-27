# (c) Lucien Mallett 2024
# This code plots the temperature response functions for AIA and saves them.

import scipy as sp
import matplotlib
import matplotlib.pyplot as plt

matplotlib.use('TkAgg')
matplotlib.rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
matplotlib.rc('text', usetex=True)

##############################################
##############################################
# INITIALIZATION
##############################################
##############################################

# Path to response functions
path = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//code//python code//"
pathresp_2011 = path + "aia_resp_110809.sav"
pathresp_2014 = path + "aia_resp_140910.sav"

savename = "C://Users//Lucien//Documents//School//Research//2023 - DEM Inversion//movies//responses.png"

##############################################
##############################################

allwaves  = ["$94 \mathrm{\AA}$",
             "$131 \mathrm{\AA}$",
             "$171 \mathrm{\AA}$",
             "$193 \mathrm{\AA}$",
             "$211 \mathrm{\AA}$",
             "$335 \mathrm{\AA}$"]
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
axs.set_ylabel(r"AIA Temp. Response $(\mathrm{DN}\;\mathrm{cm}^5 \; \mathrm{px}^{-1} \; \mathrm{s}^{-1})$")

##############################################
# 2011 response functions

# load in response functions as nx1 arrays
respinfo = sp.io.readsav(pathresp_2011, python_dict = False, verbose = False)
logte= (respinfo["aia_tresp"][0])[0]
r94  = (respinfo["aia_tresp"][0])[2]
r131 = (respinfo["aia_tresp"][0])[3]
r171 = (respinfo["aia_tresp"][0])[4]
r193 = (respinfo["aia_tresp"][0])[5]
r211 = (respinfo["aia_tresp"][0])[6]
r335 = (respinfo["aia_tresp"][0])[7]
resp = [r94, r131, r171, r193, r211, r335]

for j in range(6):
    respj = resp[j]
    axs.semilogy(logte, respj, color = linecolors[j], label = allwaves[j])

print(logte)

##############################################
# 2014 response functions

# load in response functions as nx1 arrays
respinfo = sp.io.readsav(pathresp_2014, python_dict = False, verbose = True)
logte= (respinfo["aia_tresp"][0])[0]
r94  = (respinfo["aia_tresp"][0])[2]
r131 = (respinfo["aia_tresp"][0])[3]
r171 = (respinfo["aia_tresp"][0])[4]
r193 = (respinfo["aia_tresp"][0])[5]
r211 = (respinfo["aia_tresp"][0])[6]
r335 = (respinfo["aia_tresp"][0])[7]
resp = [r94, r131, r171, r193, r211, r335]

for j in range(6):
    respj = resp[j]
    axs.semilogy(logte, respj, color = linecolors[j], label = allwaves[j], linestyle = 'dashed')

print(respinfo["aia_tresp"])

##############################################
# Plot

axs.legend(ncols = 2)
plt.tight_layout()
plt.savefig(fname = savename)
plt.show()