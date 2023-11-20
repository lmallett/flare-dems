import datetime
import matplotlib.pyplot as plt

from astropy.visualization import time_support

from sunpy import timeseries as ts
from sunpy.net import Fido
from sunpy.net import attrs as a

tstart = "2011-08-09 07:00"
tend = "2011-08-09 10:00"


# tstart = "2014-09-10 17:00"
# tend = "2014-09-10 23:00"


result = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"))
print(result)


result_goes15 = Fido.search(a.Time(tstart, tend), a.Instrument("XRS"), a.goes.SatelliteNumber(15), a.Resolution("flx1s"))
print(result_goes15)

file_goes15 = Fido.fetch(result_goes15)
goes_15 = ts.TimeSeries(file_goes15)

goes_15 = goes_15.truncate(tstart, tend)

goes_15.peek()
