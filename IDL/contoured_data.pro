;+
; Produces a series of images in outdir that display the datacube in a specific channel.
; Contoured on top of it is the saturation (red) in that channel, and the emission cube failures (black).
;
; WRITTEN BY:
;   Lucien Mallett
;
; ASSUMPTIONS:
;   dir will contain a datacube and emcube with all channels.
;   outdir
;       | 94_imgs
;       | 131_imgs
;       | 171_imgs
;       | 193_imgs
;       | 211_imgs
;       | 335_imgs
;
; INPUT:
;   (OPTIONAL) without: can be one of (94, 131, 171, 211, 335) and will produce
;       an image of data from that wavelength; otherwise it will do all of them.
;
; OUTPUT:
;   A series of images in outdir, displaying the datacube in a specific channel.
;   Contoured on top of it is the saturation (red) in that channel, and the emission cube failures (black).
;
; ADDITIONAL DOCUMENTATION:
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/image/wr_movie.pro
;   https://heliodocs.com/xdoc/xdoc_print.php?file=$SSW/gen/idl/graphics/saveimage.pro
;-

PRO contoured_data, band = band, includeem = includeem

    ; ====================
    ; INITIALIZATION
    ; ====================

    ; ; August 2011 flare
    ; xregion = [20,250]
    ; yregion = [100,400]
    ; dir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    ; IF NOT(KEYWORD_SET(includeem)) THEN BEGIN
    ;     includeem = 0
    ;     outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_data_contour_saturation/'
    ; ENDIF ELSE BEGIN
    ;     outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/imgs/imgs_data_contour_sat-em/'
    ; ENDELSE

    ; September 2014 flare
    xregion = [0,500]
    yregion = [0,400]
    dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all/'
    IF NOT(KEYWORD_SET(includeem)) THEN BEGIN
        includeem = 0
        outdir  = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_data_contour_saturation/'
    ENDIF ELSE BEGIN
        outdir  = '/home/reu23-student/Desktop/Research/2014_09_10/imgs/imgs_data_contour_sat-em/'
    ENDELSE

    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1
    
    IF NOT(KEYWORD_SET(band)) THEN band = [94, 131, 171, 193, 211, 335]

    SPAWN, 'ls ' + dir + 'aia_emcube_*.sav', filelist

    ; ====================
    ; DATA PROCESSING
    ; ====================

    ; get all data and store it into a lists of arrays

    imgs94  = LIST()
    imgs131 = LIST()
    imgs171 = LIST()
    imgs193 = LIST()
    imgs211 = LIST()
    imgs335 = LIST()

    oversat94   = LIST()
    oversat131  = LIST()
    oversat171  = LIST()
    oversat193  = LIST()
    oversat211  = LIST()
    oversat335  = LIST() 
    alloversat  = LIST()

    emcubelist  = LIST()

    FOREACH file, filelist DO BEGIN $
        RESTORE, file
        print, "Loading in ", file.substring(-27,-5)

        ; adds datacubes to list of frames for the underlying image
        cube    = TEMPORARY(datacube[xstart:xend, ystart:yend, *])

        imgs94.Add,    cube[*,*,0]
        imgs131.Add,   cube[*,*,1]
        imgs171.Add,   cube[*,*,2]
        imgs193.Add,   cube[*,*,3]
        imgs211.Add,   cube[*,*,4]
        imgs335.Add,   cube[*,*,5]

        ; adds oversaturation list to list of frames to be contoured
        oversat = TEMPORARY(satmap[xstart:xend, ystart:yend, *])

        oversat94.Add,  oversat[*,*,0]
        oversat131.Add, oversat[*,*,1]
        oversat171.Add, oversat[*,*,2]
        oversat193.Add, oversat[*,*,3]
        oversat211.Add, oversat[*,*,4]
        oversat335.Add, oversat[*,*,5]

        oversat = TOTAL(oversat,3, /INTEGER)
        alloversat.Add, oversat

        ; list where emission cubes converge / do not converge
        ; where emcubes converge, set the pixel to -1
        emcube  = TEMPORARY(emcube[xstart:xend, ystart:yend, *])
        emcubetot   = TEMPORARY(TOTAL(emcube,3))
        emcubetot[WHERE(emcubetot GT 0, /NULL)] = -1

        emcubelist.Add, emcubetot

    ENDFOREACH

    ; wavelist: list of datacube lists [{94 imgs},{131 imgs},{171 imgs},{211 imgs},{335 imgs}]
    ; measure: some sublist of [94, 131, 171, 193, 211, 335] depending on which band/s are specified
    ; oversatlist: list of satmap lists [{94 imgs},{131 imgs},{171 imgs},{211 imgs},{335 imgs}]
    wavelist = LIST()
    oversatlist = LIST()
    measure = LIST()

    ; covers case where we are asked for multiple bands
    FOR i=0, N_ELEMENTS(band)-1 DO BEGIN  $
        CASE band[i] OF
            94: BEGIN
                wavelist.Add, imgs94
                measure.Add, 94
                oversatlist.Add, oversat94
                END
            131: BEGIN
                wavelist.Add, imgs131
                measure.Add, 131
                oversatlist.Add, oversat131
                END
            171: BEGIN
                wavelist.Add, imgs171
                measure.Add, 171
                oversatlist.Add, oversat171
                END
            193: BEGIN
                wavelist.Add, imgs193
                measure.Add, 193
                oversatlist.Add, oversat193
                END
            211: BEGIN
                wavelist.Add, imgs211
                measure.Add, 211
                oversatlist.Add, oversat211
                END
            335: BEGIN
                wavelist.Add, imgs335
                measure.Add, 335
                oversatlist.Add, oversat335
                END
        ENDCASE
    ENDFOR

    ; kernel = GAUSSIAN_FUNCTION([1,1], WIDTH = 3)
    ; kernel = [[0, 1/8, 0], [1/8,1/2,1/8],[0,1/8,0]]

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
    ; plot datacubes, overlay emission cube contours

    FOR i=0, N_ELEMENTS(wavelist)-1 DO BEGIN $               ; for each wavelength

        ; load the appropriate cmap
        aia_lct, rr, gg, bb, wave = measure[i], /load
        LINECOLOR, 253, 'black'
        LINECOLOR, 254, 'red'

        waveimgs = wavelist[i]          ; data in wavelength i
        oversatimgs = oversatlist[i]    ; satmap in wavelength i

        print, 'Currently loading' + string(measure[i]) + ' A'

        FOR j=0, N_ELEMENTS(waveimgs)-1 DO BEGIN $           ; for each image in a certain wavelength, i.e., for each timestep

            data = BYTSCL((waveimgs[j])[*,*]^0.25, min = 0, TOP = 252)
            SZ = SIZE(data)

            WINDOW, 30, XSIZE = SZ[1], YSIZE = SZ[2];, /PIXMAP

            ; i don't understand how IDL handles graphic devices and
            ; at this point i'm too afraid to ask
            CONTOUR, oversatimgs[j], XSTYLE = 4, YSTYLE = 4, /NODATA, /DEVICE
            PX = !X.WINDOW; * !D.X_VSIZE
            PY = !Y.WINDOW; * !D.Y_VSIZE
            TV, data

            IF includeem THEN BEGIN
                emimg = emcubelist[j]
                CONTOUR, emimg, XSTYLE = 5, YSTYLE = 5,$
                    LEVELS = [0], $
                    POSITION = [PX[0], PY[0], PX[0]+SZ[1]-1, PY[0]+SZ[2]-1], $
                    C_COLORS = [253], $
                    /NOERASE, /DEVICE
            ENDIF

            CONTOUR, oversatimgs[j], XSTYLE = 5, YSTYLE = 5,$   ; suppress axis range + force axis range
                LEVELS = [1,3,5], $
                POSITION = [PX[0], PY[0], PX[0]+SZ[1]-1, PY[0]+SZ[2]-1], $
                C_COLORS = [254], $
                /NOERASE, /DEVICE

            trimg = tvrd(/true)
            ; wdel, 30

            filename = outdir + STRTRIM(measure[i], 2) $ ; finds folder
                + '_imgs/aia_datacube_' + STRTRIM(times[j], 2) + '_' + STRTRIM(j, 2) + '.png'
            WRITE_PNG, filename, trimg

        ENDFOR

    ENDFOR

END
