PRO process_110809_x6


;xrt files
files=find_file('/Users/awinebar/work/projects/database/110809_x6/xrt/*fits')
nfiles=n_elements(files)

read_xrt,files,ind,/nodata

m=where(ind.chip_sum NE 8 and ind.ycen GT 0 and ind.xcen GT 800 and ind.ec_fw2_ NE 'Gband')

read_xrt,files[m],i2,d2

xrt_prep,i2,d2,xindex,xdata,/normalize,/float,grade_map=grade_map,uncert_map=uncert_map

save,xindex,xdata,file='/Users/awinebar/work/projects/database/110809_x6/xrt/xrt_data.sav'

end

