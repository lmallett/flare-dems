x; This procedure takes in emission cube files in dir
; and produces emission cubes WITHOUT one channel.
; Does so for each channel and saves them as
;       aia_m@@_emcube_##-##-##T##_##_##.###.sav
; where @ is the channel (94, 131, etc.) and # is the UTC timestamp. For example,
;       outidir/m94/aia_m94_emcube_2011-08-09T07_48_15.57Z.sav
;
; WRITTEN BY:
;   Amy Winebarger
;   Edited/formatted by Lucien Mallett
;
; REQUIREMENTS:
;       M. Cheung's aia_sparse_em_solve program.
;
; ASSUMPTIONS:
;
;       dir
;           | all
;       dir/all should contain .sav files made from all channels.
;
;       outdir
;           | m94
;           | m131
;           | m171
;           | m193
;           | m211
;           | m335
;
; OUTPUT:
;       A set of .sav files in folders in outdir that contains emission cubes that were produced without one channel.
;
; ==========================================================================================

PRO without_channels

    ; dir	= '/home/reu23-student/Desktop/Research/110809_x6/emcubes/all/'
    ; outdir = '/home/reu23-student/Desktop/Research/110809_x6/emcubes/'
    ; outdir = 'E:/emcubes_110809/emcubes_m94/'
    dir = '/data/lucien/emcubes/'
    outdir = '/data/lucien/emcubes_mchannels/'

    xregion = [120,200]
    yregion = [100,250]

    time = '2011-08-09 07:48:57'

    RESTORE, 'aia_resp_110809.sav'

    nlgT = N_ELEMENTS(logt)
    lgTmin = MIN(logt)
    dlgT = ROUND((logt[1]-logt[0])*10.)/10.

    PRINT, 'Making temperature responses'
    ; creates an array, where rows represent bands and cols are logT
    tresp	= DBLARR(6,nlgt)
    tresp[0,*] = aia_tresp.resp_94
    tresp[1,*] = aia_tresp.resp_131
    tresp[2,*] = aia_tresp.resp_171
    tresp[3,*] = aia_tresp.resp_193
    tresp[4,*] = aia_tresp.resp_211
    tresp[5,*] = aia_tresp.resp_335

    ; load files
    files = FILE_SEARCH(dir + 'aia_emcube_*.sav')
    nfiles= N_ELEMENTS(files)

    wave	= [94,131,171,193,211,335]
    nwave	= N_ELEMENTS(wave)
    ww 		= INDGEN(nwave)

    ; for each file
    FOR i=0, nfiles-1 DO BEGIN

        RESTORE, files[i]
        PRINT, files[i]

        FOR j=0, nwave-1 DO BEGIN       ; iterate through each wavelength

            ; m finds the indices of the wavelengths that we want to keep
            m = WHERE(ww NE j)

            ; loads the columns of tresp that we want to keep
            trespin = tresp[m,*]
            wavein	= wave[m]
            exptime0 = exptimeout[m]

            exptimestr	= '[' + STRJOIN(STRING(exptime0, format='(F8.5)'),',') + ']'
            wavestr	    = '[' + STRJOIN(trim(wavein),',') + ']'

            tolfunc =  'aia_bp_estimate_error(y*' + exptimestr + $
                    ',[94,131,171,193,211,335],num_images=('+ strtrim(string(1^2, format = '(I4)'),2) + ')/' + exptimestr + ')'

            ; loads in data for the wavelengths that are not the wavelength we want removed
            datacubein	= datacube[*,*,m]

            ; solves the emcubes
            aia_sparse_em_init, timedepend_date = time,		$
                    use_lgtaxis = FINDGEN(nlgT)*dlgT+lgTmin,$
                    use_tresp	= trespin

            aia_sparse_em_solve, datacubein, tolfunc = tolfunc, tolfac = 1.4, $
                    oem	= emcube, 							$
                    status	= statuscube, 					$
                    coeff	= coeff,						$
                    adaptive_tolfac = adaptive_tolfac_aia,	$
                    tolmap	= tolmap_aia

            ; save
            r0 = STRPOS(files[i],'2011')
            fn = 'aia_m'+trim(wave[j])+'_emcube_'+STRMID(files[i],r0)

            print, outdir+fn
            print, "Saving ", fn
            SAVE, emcube, statuscube, logt, datacube, waveout, exptimeout, filename = outdir + fn
        ENDFOR

    ENDFOR

END


