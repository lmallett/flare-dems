;+
; Plots the number of pixels that are failing, over time.
;
; WRITTEN BY:
;   Lucien Mallett
;-


PRO numfails

    ; ====================
    ; INITIALIZATION
    ; ====================

    ; August 2011 flare
    ; xregion = [20,250]
    ; yregion = [100,400]
    xregion = [70,200]
    yregion = [200,315]
    dir     = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/other/'

    ; September 2014 flare
    ; xregion = [100,500]
    ; yregion = [50,300]
    xregion = [200,400]
    yregion = [120,200]
    dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all/'
    outdir = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/other/'

    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1

    wave = [94, 131, 171, 193, 211, 335]

    ; ====================
    ; DATA PROCESSING
    ; ====================

    tref    = '2023-01-01T00:00:00.00Z'
    sref    = anytim2tai(tref)

    ; creates a list of string arrays, one for each emcube we're dealing with
    ; each string array is the filenames for each frame 
    emfiles = FILE_SEARCH(dir + '*.sav')

    ; fills an array of the form
    ;   time    |   # black pixels
    ;=============================
    ;   frame1  |   frame1data
    ;   frame2  |   frame2data
    ;   etc     |   etc

    numdata = MAKE_ARRAY(2, N_ELEMENTS(emfiles), /LONG)

    FOR i=0, N_ELEMENTS(emfiles)-1 DO BEGIN $

        ; count failures

        print, "Loading ", (emfiles[i]).substring(-27)
        RESTORE, emfiles[i]
        emcube = TOTAL(emcube,3)^0.25
        emcube = emcube[xstart:xend, ystart:yend]

        wherepix = WHERE(emcube EQ 0, failct)

        ; timestamps

        timestamp = emfiles[i]
        timestamp = timestamp.substring(-27,-5)
        timestamp = repstr(timestamp, '_', ':')
        timestamp = anytim2tai(timestamp)

        numdata[0,i] = timestamp
        numdata[1,i] = failct

        ; update tref so that it's the minimum time we've found
        IF timestamp LE sref THEN sref = timestamp
    ENDFOR

    loadct, 8
    window, 0, xsize = xend-xstart-1, ysize = yend-ystart-1
    tv, bytscl(emcube)

    ; edits the data so that it starts at an appropriate timestamp
    tadjust = MAKE_ARRAY(2, N_ELEMENTS(emfiles), VALUE = -sref)
    tadjust[1,*]= 0
    numdata = numdata + tadjust

    mintime = min(numdata[0,*])
    maxtime = max(numdata[0,*])
    print, mintime, maxtime
    ; ====================
    ; IMAGE PROCESSING
    ; ====================
    ; colors  = ['b','g','r','c','m','hot_pink', 'b']

    p = PLOT(numdata[0,*], numdata[1,*], $
        xtitle = 'time (s)',            $
        ytitle = '# emcube failures',   $
        title  = 'Number of DEM failures', $
        xrange = [mintime, maxtime])

    p2= PLOT(indgen(N_ELEMENTS(numdata[0,*])), numdata[1,*], $
        xtitle = 'Frame #',            $
        ytitle = '# emcube failures',   $
        title  = 'Number of DEM failures',$
        xrange = [0, N_ELEMENTS(numdata[0,*])])



    p.Save, outdir + 'flare_emcube_fails.png'


END