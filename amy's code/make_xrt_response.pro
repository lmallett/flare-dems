; IDL procedure by Athiray
; Copyright (c) 2019, MaGIXS Sounding rocket mission, NASA Marshall Space Flight Center,  All rights reserved.
;       Unauthorized reproduction is allowed.


; Start		: 04 May 2023 12:37
; Last Mod 	: 06 May 2023 23:13

;-------------------  Details of the program --------------------------;
;writes out only Bethin (currently).
;time='20210730'
;xrt_tresp=make_xrt_response(time,/feldman)


function make_xrt_response, time, feldman=feldman,schmelz=schmelz,scott=scott,filter=filter
h=6.634d-34
c=2.9979d8
E_conv=(h*c/1.6d-19)*1d7

dir='xrt_responses/'


if keyword_set(feldman) then begin
    abund_name='$SSW/packages/chianti/dbase/abundance/sun_coronal_1992_feldman_ext.abund'
    outname='feldman_m'
    dbase_fname='xrt_response_database_'+outname+'.genx'    
endif

if keyword_set(schmelz) then begin
    abund_name='$SSW/packages/chianti/dbase/abundance/sun_coronal_2012_schmelz_ext.abund'
    outname='schmelz_m'
    dbase_fname='xrt_response_database_'+outname+'.genx'        
endif

if keyword_set(scott) then begin
    abund_name='$SSW/packages/chianti/dbase/abundance/sun_photospheric_2015_scott.abund'
    outname='scott_m'
    dbase_fname='xrt_response_database_'+outname+'.genx'        
endif
ioneq_name='$SSW/packages/chianti/dbase/ioneq/chianti.ioneq'
pressure=1d15

verstr=ch_get_version()

comment='Hinode/XRT responses calculated in the 1-60 Angstrom range with CHIANTI version '+$
   verstr+' data '

if n_elements(density) gt 0 then $
comment=[comment, ' at a fixed density= '+ trim(density)] else $
if n_elements(pressure) gt 0 then $
comment=[comment, ' at a fixed pressure= '+ trim(pressure)]

logt=indgen(20)*0.1+5.6
temp = 10.^logt
wmin=1.
wmax=60.

if not file_exist(dir+dbase_fname) then begin
    isothermal, wmin, wmax, 0.1, temp, lambda,spectrum,list_wvl,list_ident,$
	pressure=pressure, edensity=density,/photons,/cont , $
	abund_name=abund_name, ioneq_name=ioneq_name, verbose=verbose
     database={abundance:outname, lambda:lambda,spectrum:spectrum,$
	list_wvl:list_wvl,list_ident:list_ident,pressure:pressure}
    savegen,database,file=dir+dbase_fname
endif


restgen,dbase,file=dir+dbase_fname
lambda=dbase.lambda
spectrum=dbase.spectrum

energy=E_conv/lambda
elecperphoton=energy/0.00365

comment=[comment, ' with the elemental abundance file: '+abund_name]
comment=[comment, ' with the ionisation equilibrium file: '+ioneq_name]
wave_resp=make_xrt_wave_resp(contam_time=time)

;Index for filters 
;bethin=4, bemed=5,bethick=8
;if filter EQ 'Be-thin' then m = 4 ELSE stop

m=where(trim(wave_resp.name) EQ trim(filter))
wave=wave_resp.effar.wave
xrt_ea=wave_resp.effar.eff_area

filter_ea=xrt_ea[*,m[0]]
xrt_wave=wave[*,0]

eff_filter=interpol(filter_ea[0:3985], xrt_wave[0:3985], lambda)>0.
sterad_xrt_pix = 1.02860*1.02860/4.25d10
electronperDN=59.0
sp_conv= spectrum & sp_conv[*,*]=0.
for i=0,n_elements(temp)-1 do $
    sp_conv[*, i]= sterad_xrt_pix*spectrum[*, i]  * eff_filter *elecperphoton/electronperDN
resp_filter=total(sp_conv,1)

tresp={temperatures:logt,comment: comment, $
    resp_filter:resp_filter,filter:filter}
;ofname='xrt_resp_'+time+'_'+outname+'.sav'
;save,tresp,file=ofname

return,tresp
END

