PRO decon_aia

common thispsf,psf


wave=['94','131','171','193','211','335']
nwave=n_elements(wave)

;the directory where the data is:
dir='/Users/awinebar/work/projects/database/110809_x6/aia/'
;the center of the flare
xx=3499
yy=2422


for i=0,nwave-1 do begin

psf=aia_calc_psf(float(wave[i]))

files=find_file(dir+'*'+wave[i]+'A*fits')

nfiles=n_elements(files)
c=0

for ii=0,nfiles-1 do begin

	read_sdo,files[ii],index0,data0
	satmap0=data0*0.
	mm=where(data0 GT 15000.,nmm)
	if nmm GT 0 THEN satmap0[mm]=1
	data0 = aia_deconvolve_richardsonlucy(float(data0), psf)
	exptime0=index0.exptime
	aia_prep,index0,data0,index1,data1,/norm
	index=index1
	data=data1[xx-500:xx+500,yy-500:yy+500]
	satmap=satmap0[xx-500:xx+500,yy-500:yy+500]
	exptime=exptime0
	r1=strpos(files[ii],'.fits')
	r0=strpos(files[ii],'/',/reverse_search)
	newfilename=strmid(files[ii],0,r1)+'_decon_norm.sav'
	save,index,data,exptime,satmap,file=newfilename

endfor

endfor
end



