PRO drive_second_em_calculation

time='2011-08-09 07:48:57'
aia_tresp=make_aia_response(time,/feldman)
xrt_tresp=make_xrt_response(time,/feldman,filter='Ti-poly')

logt = aia_tresp.temperatures
nlgT = n_elements(logt)
lgTmin = min(logt)
dlgT = round((logt[1]-logt[0])*10.)/10.

        tresp=dblarr(7,nlgt)
        tresp[0,*] = aia_tresp.resp_94
	tresp[1,*] = aia_tresp.resp_131
	tresp[2,*] = aia_tresp.resp_171
	tresp[3,*] = aia_tresp.resp_193
	tresp[4,*] = aia_tresp.resp_211
	tresp[5,*] = aia_tresp.resp_335
	tresp[6,*] = xrt_tresp.resp_filter


;firs
files = find_file('/Users/awinebar/work/projects/database/110809_x6/emcubes/aia_emcube_*.sav') 
nfiles= n_elements(files)

thresh = 40

wave = [94,131,171,193,211,335]
nwave= n_elements(wave)
ww = indgen(nwave)

for i=0, nfiles-1 do begin

restore, files[i]

tem = total(emcube,3)
m = where(tem[120:200,120:250] EQ 0, nm)

if nm GT thresh then begin

	for j=0,nwave-1 do begin

	m=where(ww NE j)

	trespin = tresp[m,*]
	wavein = wave[m]
	exp = exptimeout[m]

	exptimestr = '[' + strjoin(string(exp,format = '(F8.5)'),',') + ']'
	wavestr= '[' + strjoin(trim(wavein),',') + ']'

        tolfunc = 'aia_bp_estimate_error(y*' + exptimestr + ', ' + wavestr + ')/' + exptimestr

	datacubein = datacube[*,*,m]

	aia_sparse_em_init, timedepend_date = time,$
    		use_lgtaxis = findgen(nlgT)*dlgT+lgTmin,$
    		use_tresp = trespin

	aia_sparse_em_solve, datacubein, tolfunc = tolfunc, tolfac=1.4, $
    		oem = emcube, status = statuscube, coeff = coeff , adaptive_tolfac = adaptive_tolfac_aia,$
	    	tolmap = tolmap_aia

	outdir = '/Users/awinebar/work/projects/database/110809_x6/emcubes/'
	r0 = strpos(files[i],'2011')
	fn = 'aia_m'+trim(wave[j])+'_emcube_'+strmid(files[i], r0)

	SAVE, emcube, statuscube, logt,datacube,waveout,exptimeout,FILENAME=outdir+fn
	endfor
endif

endfor

end


