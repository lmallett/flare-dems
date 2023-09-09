;+
; Produces a series of images in outdir, displaying the data for a channel at each pixel.
;
; WRITTEN BY:
;   Lucien Mallett
;
; ASSUMPTIONS:
;   dir will contain .sav files with a datacube with all channels.
;   outdir
;       | oversat_94
;       | oversat_131
;       | oversat_171
;       | oversat_193
;       | oversat_211
;       | oversat_335
;
; INPUT:
;   (OPTIONAL) without: can be one of (94, 131, 171, 211, 335) and will produce
;       an image of data from that wavelength; otherwise it will do all of them.
;
; OUTPUT:
;   A series of images in outdir, displaying the data at each pixel.
;
; ADDITIONAL DOCUMENTATION:
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/image/wr_movie.pro
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/graphics/saveimage.pro
;-

PRO imgs_datacubes, band = band

    ; ====================
    ; INITIALIZATION
    ; ====================

    ; ; August 2011 flare
    ; xregion = [20,250]
    ; yregion = [100,400]
    ; dir     = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    ; outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_datacubes/'

    ; September 2014 flare
    xregion = [0,500]
    yregion = [0,400]
    dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all/'
    outdir = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_datacubes/'

    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1

    IF NOT(KEYWORD_SET(band)) THEN band = [94, 131, 171, 193, 211, 335]

    ; ====================
    ; DATA PROCESSING
    ; ====================

    ; get all data and store it into a list of arrays (# images, x, y)

    w94 = LIST()
    w131= LIST()
    w171= LIST()
    w193= LIST()
    w211= LIST()
    w335= LIST()
    
    SPAWN, 'ls ' + dir + 'aia_emcube_*.sav', filelist

    FOREACH file, filelist DO BEGIN $
        RESTORE, file
        print, file.substring(-27,-5)

        cube = datacube[xstart:xend, ystart:yend, *]^0.25
        w94.Add,    cube[*,*,0]
        w131.Add,   cube[*,*,1]
        w171.Add,   cube[*,*,2]
        w193.Add,   cube[*,*,3]
        w211.Add,   cube[*,*,4]
        w335.Add,   cube[*,*,5]
    ENDFOREACH
    
    ; covers case where we are asked for multiple bands
    imglist = LIST()
    measure = LIST()
    FOR i=0, N_ELEMENTS(band)-1 DO BEGIN  $
        CASE band[i] OF
            94: BEGIN
                imglist.Add, w94
                measure.Add, 94
                END
            131: BEGIN
                imglist.Add, w131
                measure.Add, 131
                END
            171: BEGIN
                imglist.Add, w171
                measure.Add, 171
                END
            193: BEGIN
                imglist.Add, w193
                measure.Add, 193
                END
            211: BEGIN
                imglist.Add, w211
                measure.Add, 211
                END
            335: BEGIN
                imglist.Add, w335
                measure.Add, 335
                END
        ENDCASE
    ENDFOR

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

    FOR i=0, N_ELEMENTS(imglist)-1 DO BEGIN $               ; for each wavelength
        
        ; load the appropriate cmap
        aia_lct, rr, gg, bb, wave = measure[i], /load
        ctable  = [[rr],[gg],[bb]]

        print, 'Currently loading' + string(measure[i]) + ' A'
        waveimg = imglist[i]

        FOR j=0, N_ELEMENTS(waveimg)-1 DO BEGIN $           ; for each image in the wavelength

            data = BYTSCL(waveimg[j,*,*])
            window, 30, xsize = (xend-xstart+1), ysize = (yend-ystart+1), /pixmap
            tv, data

            trimg = tvrd(/true)
            wdel, 30

            filename = outdir + STRTRIM(measure[i], 2) $    ; finds folder
                + '_imgs/aia_datacube_' + STRTRIM(times[j], 2) + '_' + STRTRIM(j, 2) + '.png'

            WRITE_PNG, filename, trimg

        ENDFOR

    ENDFOR

END