; IDL procedure by Athiray
; Copyright (c) 2019, MaGIXS Sounding rocket mission, NASA Marshall Space Flight Center,  All rights reserved.
;       Unauthorized reproduction is allowed.

; Start		: 04 May 2023 10:39
; Last Mod 	: 06 May 2023 12:47

;-------------------  Details of the program --------------------------;
; REQUIREMENTS: 
;       CHIANTI package
;       sdo/aia/idl/util/aia_sparse_em_init.pro
;
; ASSUMPTIONS:
;   The directory with this file should contain a folder aia_responses/.
;
; INPUTS:
;       NONE
;
; OUTPUTS:
;       A .genx file in the aia response folder.
;       Returns a dictionary, tresp.
;


; AIA_tresp = make_aia_response(20210730,/feldman)

FUNCTION make_aia_response, time, feldman=feldman,schmelz=schmelz,scott=scott

    dir = 'aia_responses/'

    IF KEYWORD_SET(feldman) THEN BEGIN
        abund_name  ='$SSW/packages/chianti/dbase/abundance/sun_coronal_1992_feldman_ext.abund'
        outname ='feldman_m'
        dbase_fname ='aia_response_database_' + outname + '.genx'
    ENDIF

    IF KEYWORD_SET(schmelz) THEN BEGIN
        abund_name  ='$SSW/packages/chianti/dbase/abundance/sun_coronal_2012_schmelz_ext.abund'
        outname ='schmelz_m'
        dbase_fname ='aia_response_database_' + outname + '.genx'    
    ENDIF

    IF KEYWORD_SET(scott) THEN BEGIN
        abund_name  ='$SSW/packages/chianti/dbase/abundance/sun_photospheric_2015_scott.abund'
        outname ='scott_m'
        dbase_fname ='aia_response_database_' + outname + '.genx'    
    ENDIF
    ioneq_name  = '$SSW/packages/chianti/dbase/ioneq/chianti.ioneq'
    pressure    = 1d15

    ssw_path, /aia

    ; Solid angle of one AIA pixel:
    sterad_aia_pix = 8.46158d-12

    ; Get time-dependent, eve corrected, AIA effective areas
    aia_resp = aia_get_response(/dn, /area, /eve, timedepend_date=time)

    verstr  = ch_get_version()

    comment = 'AIA responses calculated in the 25-900 Angstrom range with CHIANTI version '+$
    verstr + ' data '


    IF N_ELEMENTS(density) GT 0 THEN $
    comment = [comment, ' at a fixed density= ' + trim(density)] ELSE $
    IF N_ELEMENTS(pressure) GT 0 THEN $
    comment = [comment, ' at a fixed pressure= '+ trim(pressure)]

    logt=indgen(20)*0.1+5.6
    temp = 10.^logt
    wmin = 25.
    wmax = 900
    resp_94 = dblarr(N_ELEMENTS(logt))
    resp_131= dblarr(N_ELEMENTS(logt))
    resp_171= dblarr(N_ELEMENTS(logt))
    resp_193= dblarr(N_ELEMENTS(logt))
    resp_211= dblarr(N_ELEMENTS(logt))
    resp_335= dblarr(N_ELEMENTS(logt))

    IF NOT file_exist(dir+dbase_fname) THEN BEGIN
        isothermal, wmin, wmax, 0.1, temp, lambda, spectrum, list_wvl, list_ident,$
        pressure = pressure, $
        edensity=density,/photons,/cont , $
        abund_name = abund_name, $
        ioneq_name=ioneq_name, verbose=verbose

        database={ $
            abundance:outname,  $
            lambda:lambda,  $
            spectrum:spectrum,  $
            list_wvl:list_wvl,  $
            list_ident:list_ident,  $
            pressure:pressure}

        savegen, database,file=dir+dbase_fname
    ENDIF

    restgen, dbase, file = dir+dbase_fname
    lambda  = dbase.lambda
    spectrum= dbase.spectrum

    comment = [comment, ' with the elemental abundance file: '+abund_name]
    comment = [comment, ' with the ionisation equilibrium file: '+ioneq_name]

    eff_94  = interpol(aia_resp.a94.ea, aia_resp.a94.wave,  lambda)
    eff_131 = interpol(aia_resp.a131.ea, aia_resp.a131.wave,  lambda)
    eff_171 = interpol(aia_resp.a171.ea, aia_resp.a171.wave,  lambda)
    eff_193 = interpol(aia_resp.a193.ea, aia_resp.a193.wave,  lambda)
    eff_211 = interpol(aia_resp.a211.ea, aia_resp.a211.wave,  lambda)
    eff_335 = interpol(aia_resp.a335.ea, aia_resp.a335.wave,  lambda)

    sp_conv = spectrum & sp_conv[*,*] = 0.
    FOR i=0,N_ELEMENTS(temp)-1 DO $
        sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_94
    resp_94 = TOTAL(sp_conv,1)

    sp_conv = spectrum & sp_conv[*,*] = 0.
    FOR i=0,N_ELEMENTS(temp)-1 DO $
        sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_131
    resp_131= TOTAL(sp_conv,1)

    sp_conv = spectrum & sp_conv[*,*] = 0.
    FOR i=0,N_ELEMENTS(temp)-1 DO $
        sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_171
    resp_171= TOTAL(sp_conv,1)

    sp_conv = spectrum & sp_conv[*,*] = 0.
    FOR i=0,N_ELEMENTS(temp)-1 DO $
        sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_193
    resp_193= TOTAL(sp_conv,1)

    sp_conv = spectrum & sp_conv[*,*] = 0.
    FOR i=0,N_ELEMENTS(temp)-1 DO $
        sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_211
    resp_211= TOTAL(sp_conv,1)

    sp_conv = spectrum & sp_conv[*,*]=0.
    FOR i=0,N_ELEMENTS(temp)-1 DO $
        sp_conv[*, i]= sterad_aia_pix*spectrum[*, i] * eff_335
    resp_335= TOTAL(sp_conv,1)

    tresp   = {temperatures:logt,comment: comment, $
        resp_94:resp_94,resp_131:resp_131,resp_171:resp_171, $
        resp_193:resp_193,resp_211:resp_211, resp_335:resp_335}

    RETURN, tresp
END