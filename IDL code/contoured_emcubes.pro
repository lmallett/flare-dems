;+
; Produces a series of images in outdir, displaying the sum total of the DEM at each pixel.
; Overlaid will be a contour displaying where data is saturated in 1, 3, and 5 channels.
;
; WRITTEN BY:
;   Lucien Mallett
;
; ASSUMPTIONS:
;   dir will contain .sav files with emcubes produced from all channels.
;
; OUTPUT:
;   A series of images in outdir, displaying the sum total of the DEM at each pixel,
;   with a contour displaying where data is saturated.
;
;   oversat array:
;       0: no channel saturated
;       1: 1  channel saturated
;       2: 2 channels saturated
;       3: 3 channels saturated
;       4: 4 channels saturated
;       5: 5 channels saturated
;       6: 6 channels saturated
;-

PRO contoured_emcubes

    ; ====================
    ; INITIALIZATION
    ; ====================

    ; ; August 2011 flare
    ; xregion = [20,250]
    ; yregion = [100,400]
    ; dir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    ; outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_contouredsat/'

    ; September 2014 flare
    xregion = [0,500]
    yregion = [0,400]
    dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all/'
    outdir  = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_contouredsat/'

    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1

    SPAWN, 'ls ' + dir + 'aia_emcube_*.sav', filelist

    ; ====================
    ; DATA PROCESSING
    ; ====================

    oversatlist = LIST()
    emcubelist  = LIST()

    FOREACH file, filelist DO BEGIN $
        RESTORE, file
        print, "Loading in ", file.substring(-27,-5)

        emcube      = TEMPORARY(emcube[xstart:xend, ystart:yend, *])
        emcubetot   = TEMPORARY(TOTAL(emcube,3))
        emcubelist.Add, emcubetot

        ; produce a 2D array with values ranging from 0 to 6
        ; where the number tells us how many channels were oversaturated in the data
        oversat     = TEMPORARY(satmap[xstart:xend, ystart:yend, *])
        oversat = TOTAL(oversat,3, /INTEGER)
        oversatlist.Add, oversat

    ENDFOREACH

    y = WHERE(oversat EQ 0, county)
    z = WHERE(emcubetot EQ 0, countz)
    print, county, countz

    ; set up timestamps
    times = LIST()
    FOR i=0, N_ELEMENTS(filelist)-1 DO BEGIN    $
        timestamp   = filelist[i]
        timestamp   = timestamp.substring(-27,-5)   ; get UTC time from filename
        times.Add, timestamp
    ENDFOR

    ; ====================
    ; IMAGE PROCESSING
    ; ====================
    
    ; colors
    
    loadct, 0

    ; max linecolor = 254
    ; LINECOLOR, 252, VALUES = [237, 248, 177]
    ; LINECOLOR, 253, VALUES = [127,205,187]
    ; LINECOLOR, 254, VALUES = [44,127,184]

    LINECOLOR, 252, VALUES = [255, 20, 236]
    LINECOLOR, 253, VALUES = [143,10,209]
    LINECOLOR, 254, VALUES = [30,0,182]

    ; create images
    FOR i=0, N_ELEMENTS(oversatlist)-1 DO BEGIN $   ; for each timestep

        data = BYTSCL((emcubelist[i])[*,*]^0.25, min =0, TOP = 251)     ; top = min linecolor used - 1. eg linecolor 253, so top is 252

        SZ = SIZE(data)

        WINDOW, 30, XSIZE = SZ[1], YSIZE = SZ[2];, /PIXMAP
        CONTOUR, oversatlist[i], XSTYLE = 4, YSTYLE = 4, /NODATA, /DEVICE

        PX = !X.WINDOW; * !D.X_VSIZE
        PY = !Y.WINDOW; * !D.Y_VSIZE
        TV, data

        CONTOUR, oversatlist[i], XSTYLE = 5, YSTYLE = 5,$
            LEVELS = [1, 3, 5], $
            POSITION = [PX[0], PY[0], PX[0]+SZ[1]-1, PY[0]+SZ[2]-1], $
            C_COLORS = [252,253,254], $
            /NOERASE, /DEVICE

        trimg = tvrd(/true)
        wdel, 30

        filename = outdir + 'aia_emcube_contouredsat_' + STRTRIM(times[i], 2) + '_' + STRTRIM(i, 2) + '.png'
        WRITE_PNG, filename, trimg

    ENDFOR

END