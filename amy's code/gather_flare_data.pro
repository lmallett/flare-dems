PRO gather_flare_data

restore,'/Users/awinebar/work/projects/database/110809_x6/xrt/xrt_data.sav'

wave=[94,131,171,193,211,335]
nwave=n_elements(wave)

i=0
files=find_file('/Users/awinebar/work/projects/database/110809_x6/aia/*'+trim(wave[i])+'A*decon_norm.sav')
nfiles=n_elements(files)

fn=strarr(nfiles,nwave)
times=strarr(nfiles,nwave)
nsat=intarr(nfiles,nwave)

for i=0,nwave-1 do begin
files=find_file('/Users/awinebar/work/projects/database/110809_x6/aia/*'+trim(wave[i])+'A*decon_norm.sav')
nfiles=n_elements(files)
	fn[*,i]=files
	for j = 0,nfiles-1 do begin
		restore,files[j]
		times[j,i]=index.t_obs
		m=where(satmap EQ 1,nm)
		nsat[j,i]=nm
	endfor
endfor

	ax=findgen(1001)*0.6
	nxx=round(1001.*0.6/xindex[0].cdelt1)
	xx=findgen(nxx)*xindex[0].cdelt1

for j=0,nfiles-1 do begin
	data0=fltarr(1001,1001,6)
	satmap0=fltarr(1001,1001,6)
	exptime0=fltarr(6)
	;read in the files associated with this time
	for i=0,nwave-1 do begin
		restore,fn[j,i]
		data0[*,*,i]=data
		satmap0[*,*,i]=satmap
		exptime0[i]=exptime
	endfor

;regrid the aia data only xrt bin size

	newdata=fltarr(nxx,nxx,nwave)
	for i=0,nwave-1 do newdata[*,*,i]=interp2d(data0[*,*,i],ax,ax,xx,xx,/grid)

	newsatmap=fltarr(nxx,nxx,nwave)
	for i=0,nwave-1 do newsatmap[*,*,i]=interp2d(satmap0[*,*,i],ax,ax,xx,xx,/grid)
	mm=where(newsatmap GT 0,nmm)
	if nmm GT 0 then newsatmap[mm]=1.

;find out if there is any xrt data +-30 s from average AIA time

	;min/max time
	dt=difftime(times[j,0],times[j,*])
	mintime=times[j,where(dt EQ min(dt))]
	maxtime=times[j,where(dt EQ max(dt))]
	dttot=max(dt)-min(dt)
	avgtime=addtime(mintime,dttot/2.)

	dtx0=difftime(addtime(avgtime,-30),xindex.date_obs)
	dtx1=difftime(addtime(avgtime,30),xindex.date_obs)
	m=where(dtx0 GE 0 and dtx1 LE 0 and xindex.ec_fw2_ EQ 'Ti_poly',nm)

	IF nm GT 0 THEN BEGIN

	if nm GT 1 THEN BEGIN
		thisxdata=xdata[*,*,m]
		thisexptime=xindex[m].exptime
		for jj=0,nm-1 do thisxdata[*,*,jj]=thisxdata[*,*,jj]*thisexptime[jj]
		thisxdata=total(thisxdata,3)/total(thisexptime)
		thisexptime=total(thisexptime)
	endif else begin
		thisxdata=xdata[*,*,m]
		thisexptime=xindex[m].exptime
	endelse


;align 94 data to xrt, aligning the quarter root gets the limb better

	nxrt=n_elements(xdata[*,0])

	diffx=0
	diffy=0
	offset= 1000

	while max(abs(offset)) GT 1 DO BEGIN

		suba=newdata[nxx/2.-nxrt/2.+diffx:nxx/2.+nxrt/2.-1+diffx,nxx/2.-nxrt/2.+diffy:nxx/2.+nxrt/2.-1+diffy,0]
		ss=fltarr(nxrt,nxrt,2)
		ss[*,*,0]=(suba>0)^.25
		ss[*,*,1]=(thisxdata>0)^.25
		offset=tr_get_disp(ss)
		diffx=diffx+offset[0,1]
		diffy=diffy+offset[1,1]
		print,max(abs(offset))

	endwhile

	datacube=fltarr(nxrt,nxrt,7)
	datacube[*,*,0:5]=newdata[nxx/2.-nxrt/2.+diffx:nxx/2.+nxrt/2.-1+diffx,nxx/2.-nxrt/2.+diffy:nxx/2.+nxrt/2.-1+diffy,0:5]>0.
	datacube[*,*,6]=thisxdata>0.
	satmap=fltarr(nxrt,nxrt,7)
	satmap[*,*,0:5]=newsatmap[nxx/2.-nxrt/2.+diffx:nxx/2.+nxrt/2.-1+diffx,nxx/2.-nxrt/2.+diffy:nxx/2.+nxrt/2.-1+diffy,0:5]
	waveout=[trim(string(wave)),'Ti_poly']
	exptimeout=[exptime0,thisexptime]


	endif ELSE BEGIN

	;if there is no xrt data, just write out the AIA data using the last XRT settings
        datacube=newdata[nxx/2.-nxrt/2.+diffx:nxx/2.+nxrt/2.-1+diffx,nxx/2.-nxrt/2.+diffy:nxx/2.+nxrt/2.-1+diffy,*]>0.
	satmap=newsatmap[nxx/2.-nxrt/2.+diffx:nxx/2.+nxrt/2.-1+diffx,nxx/2.-nxrt/2.+diffy:nxx/2.+nxrt/2.-1+diffy,*]

        waveout=[trim(string(wave))]
        exptimeout=[exptime0]
	ENDELSE

	

;make a new filename
	outdir='/Users/awinebar/work/projects/database/110809_x6/datacubes/'
	r0=strpos(fn[j,0],'2011')
	r1=strpos(fn[j,0],'.image')
	filename=strmid(fn[j,0],r0,r1-r0)
;save the data cube
	save,datacube,waveout,exptimeout,file=outdir+filename+'.sav'
endfor

end
	

