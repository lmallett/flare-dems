;+
; Produces a series of images in outdir, displaying the statuscube result at each pixel.
;
; WRITTEN BY:
;   Lucien Mallett
;
; ASSUMPTIONS:
;   dir will contain .sav files with a statuscube.
;
; INPUT:
;   (OPTIONAL) without: can be one of (94, 131, 171, 211, 335) and will produce
;       an image of the statuscube that was calculated without that wavelength.
;
; OUTPUT:
;   A series of images in outdir, displaying the status of the DEM inversion at each pixel.
;       COLOR         | ERROR CODE: Explanation
;       ==============================================================
;       Light blue:   | 1: The objective function is unbounded.
;       Dark blue:    | 2: No solution satisfies the given constraints.
;       Light green:  | 3: The routine did not converge.
;       Dark green:   | 10: SIMPLEX very occasionally returns negative coeffs even though positivity is a constraint. Now these are set to zero.
;
; ADDITIONAL DOCUMENTATION:
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/image/wr_movie.pro
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/graphics/saveimage.pro
;-

PRO imgs_statuscubes, without = band, nolegend = nolegend

    ; ====================
    ; INITIALIZATION
    ; ====================

    ; August 2011 flare
    xregion = [20,250]
    yregion = [100,400]
    dir     = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_statuscubes/'
    
    ; ; September 2014 flare
    ; xregion = [0,500]
    ; yregion = [0,400]
    ; dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all/'
    ; outdir = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_statuscubes/'

    IF NOT(KEYWORD_SET(band)) THEN BEGIN
        SPAWN, 'ls ' + dir + 'aia_emcube_*.sav', filelist
        band = 'all'
    ENDIF ELSE BEGIN
        band = 'm' + STRTRIM(band, 2)
        SPAWN, 'ls ' + dir + 'aia_'+ STRING(band) +'_emcube_*.sav', filelist
    ENDELSE

    IF NOT(KEYWORD_SET(nolegend)) THEN nolegend = 0

    ; ====================
    ; DATA PROCESSING
    ; ====================

    ; get all data and store it into a list of arrays (# images, x, y)
    statcubes = LIST()

    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1

    FOREACH file, filelist DO BEGIN $
        RESTORE, file
        statcube = statuscube[xstart:xend, ystart:yend]
        statcubes.Add, statcube
    ENDFOREACH

    ; set up timestamps
    times = LIST()
    FOR i=0, N_ELEMENTS(filelist)-1 DO BEGIN $
        timestamp   = filelist[i]
        timestamp   = timestamp.substring(-27,-5)   ; get UTC time from filename
        times.Add, timestamp
    ENDFOR

    ; ====================
    ; IMAGE PROCESSING
    ; ====================

    ; deal with colors

    numlevels = CEIL(!D.TABLE_SIZE/6)
    level = INDGEN(!D.TABLE_SIZE)/numlevels
    
    ; ctable = $
    ;     [[0,0,0], $
    ;     [255,255,204], $     
    ;     [161,218,180], $
    ;     [65,182,196], $
    ;     [34,94,168]]

    ; R = ctable[0, level]
    ; G = ctable[1, level]
    ; B = ctable[2, level]
    ; R[!D.TABLE_SIZE-1] = 34
    ; G[!D.TABLE_SIZE-1] = 94
    ; B[!D.TABLE_SIZE-1] = 168

    ctable = $
        [[0,0,0], $
        [166,206,227],$     
        [31,120,180], $
        [178,223,138],  $
        [51,160,44],$
        [235, 52, 183]]

    R = ctable[0, level]
    G = ctable[1, level]
    B = ctable[2, level]
    R[!D.TABLE_SIZE-1] = 235
    G[!D.TABLE_SIZE-1] = 52
    B[!D.TABLE_SIZE-1] = 183

    tvlct, R, G, B
    
    ; loadct, 38
    ; xlen = xend-xstart +1
    ; clen = CEIL(xlen/6)


    FOR i=0, N_ELEMENTS(statcubes)-1 DO BEGIN $

        print, times[i]
    
        img = statcubes[i]

        ; img[0:clen,0] = 0
        ; img[clen*1:clen*2,0] = 1
        ; img[clen*2:clen*3,0] = 2
        ; img[clen*3:clen*4,0] = 3
        ; img[clen*4:clen*5,0] = 10
        ; img[clen*5:*,0] = 11

        ; changes values, for color scaling
        img[WHERE(img EQ 10, /NULL)] = 4
        img[WHERE(img EQ 11, /NULL)] = 5

        window, 30, xsize = (xend-xstart+1), ysize = (yend-ystart+1), /pixmap
        tv, bytscl(img[*,*], min = 0, max = 5)

        IF NOT nolegend THEN BEGIN
            AL_LEGEND, ['0','1','2','3','10','11'], $
                /fill, $
                psym = 15 + intarr(6), $
                colors = [0,52,105,155,180,240], $
                background_color = 'white', $
                pos = [-0.1,1.05], $
                charsize = 1.25
        ENDIF

        trimg = tvrd(/true)
        wdel, 30

        filename = outdir + times[i] + STRING(i) + '.png'
        WRITE_PNG, STRCOMPRESS(filename), trimg

    ENDFOR
 
END
