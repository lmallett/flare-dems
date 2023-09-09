;+
; Produces a series of images in outdir, displaying the sum total of the DEM at each pixel.
; Produces images from emcube save files.
;
; WRITTEN BY:
;   Lucien Mallett
;
; ASSUMPTIONS:
;   Directory structure is as follows:
;   dir
;       | all
;       | m94
;       | m131
;       | m171
;       | m193
;       | m211
;       | m335
;   dir/all will contain .sav files with emcubes produced from all channels.
;
; INPUT:
;   (OPTIONAL) without: can be one of (94, 131, 171, 211, 335) and will produce
;       an image of the emcube that was calculated without that wavelength.
;
; OUTPUT:
;   A series of images in outdir, displaying the sum total of the DEM at each pixel.
;
; ADDITIONAL DOCUMENTATION:
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/image/wr_movie.pro
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/graphics/saveimage.pro
;
; There will be a memory leak if you try to do the object-oriented SAVE IMAGE method. Thanks IDL.
;-

PRO imgs_emcubes, without = band

    ; ====================
    ; INITIALIZATION
    ; ====================

    ; ;;;;;;;;
    ; ; August 2011 flare
    ; xregion = [20,250]
    ; yregion = [100,400]
    ; dir     = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/'
    ; outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_emcubes/'

    ; IF NOT(KEYWORD_SET(band)) THEN  BEGIN
    ;     loadct, 8 ; 3 1
    ;     outdir = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_emcubes/all/'
    ;     SPAWN, ('ls ' + dir + 'all/*.sav'), filelist
    ;     band = 'all'
    ; ENDIF ELSE BEGIN
    ;     aia_lct, rr, gg, bb, wave = band, /load
    ;     band = 'm' + STRTRIM(band, 2)
    ;     SPAWN, ('ls ' + dir + band + '/aia_'+ band +'_emcube_*.sav'), filelist
    ;     outdir = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_emcubes/' + band + '/'
    ; ENDELSE
    ; ;;;;;;;;

    ;;;;;;;
    ; September 2014 flare
    xregion = [0,500]
    yregion = [0,400]
    dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/'
    outdir = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_emcubes/'

    IF NOT(KEYWORD_SET(band)) THEN  BEGIN
        loadct, 8 ; 3 1
        outdir = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_emcubes/'
        SPAWN, ('ls ' + dir + 'all/*.sav'), filelist
        band = 'all'
    ENDIF ELSE BEGIN
        aia_lct, rr, gg, bb, wave = band, /load
        band = 'm' + STRTRIM(band, 2)
        SPAWN, ('ls ' + dir + band + '/aia_'+ band +'_emcube_*.sav'), filelist
        outdir = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_emcubes/' + band + '/'
    ENDELSE
    ;;;;;;;

    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1
    
    ; ====================
    ; DATA PROCESSING
    ; ====================

    print, "Restoring files..."

    ; get all data and store it into a list of arrays (# images, x, y)
    ; creates an image of the summed DEMs; appends to a list of images
    emcubes = LIST()
    FOREACH file, filelist DO BEGIN         $
        RESTORE, file
        emcube = (TOTAL(emcube,3))^0.2
        emcube = emcube[xstart:xend, ystart:yend]
        emcubes.Add, emcube
    ENDFOREACH

    ; set up timestamps
    times = LIST()
    FOR i=0, N_ELEMENTS(filelist)-1 DO BEGIN    $
        timestamp   = filelist[i]
        timestamp   = timestamp.substring(-27,-5)   ; get UTC time from filename
        times.Add, timestamp
    ENDFOR

    ; ====================
    ; PRODUCE IMAGES
    ; ====================

    print, "Making images..."
    FOR i=0, N_ELEMENTS(emcubes)-1 DO BEGIN $

        print, times[i]
        data = BYTSCL(emcubes[i,*,*], min = 0)

        window, 30, xsize = (xend-xstart+1), ysize = (yend-ystart+1), /pixmap
        tv, data
        ;xyouts, 0,0, times[i], charsize = 1.5
        ;xyouts, 0,0.95, band, charsize = 1.5, /NORMAL
        trimg = tvrd(/true)
        wdel, 30
        
        ; im = IMAGE(         $
        ;     bytscl(data),   $
        ;     DIMENSIONS = [(xend-xstart+1), (yend-ystart+1)], $
        ;     RGB_TABLE = ctable,  $
        ;     MARGIN = 0,     $
        ;     /BUFFER)

        filename = outdir + band + '_' + times[i] + STRING(i) + '.png'
        ; im.Save, STRCOMPRESS(filename)

        WRITE_PNG, STRCOMPRESS(filename), trimg

    ENDFOR

END