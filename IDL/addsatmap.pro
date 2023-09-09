PRO addsatmap

    badrawdir  = '/home/reu23-student/Desktop/Research/2011_08_09/data_aia/decon_norm_badsat/'
    goodrawdir = '/home/reu23-student/Desktop/Research/2011_08_09/data_aia/decon_norm/'

    bademdir   = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all-badmap/'
    datadir = '/home/reu23-student/Desktop/Research/2011_08_09/datacubes/'
    satdir  = '/home/reu23-student/Desktop/Research/2011_08_09/data_aia/satmaps/'

    goodemdir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'

    emfiles     = FILE_SEARCH(bademdir + '*.sav')
    datafiles   = FILE_SEARCH(datadir + '*.sav')

    wave = ['94A', '131A', '171A','193A', '211A', '335A']

    x0 = 350
    x1 = 650
    y0 = 275
    y1 = 775

    ; ====================
    ; For correcting decon_norm files.
    ; ====================

    ; FOR i=0, N_ELEMENTS(wave)-1 DO BEGIN

    ;     rawdatafiles= FILE_SEARCH(badrawdir + wave[i] + '/*.sav')
    ;     satfiles    = FILE_SEARCH(satdir +  wave[i] + '/*.sav')

    ;     FOR j=0, N_ELEMENTS(rawdatafiles)-1 DO BEGIN

    ;         RESTORE, rawdatafiles[j]
    ;         RESTORE, satfiles[j]
    ;         satmap = satmap[x0:x1, y0:y1]

    ;         r0 = STRPOS(rawdatafiles[j], '/', /reverse_search)
    ;         fn = rawdatafiles[j].substring(r0)
    ;         print, goodrawdir + wave[i] + fn
    ;         SAVE, data, exptime, index, satmap, FILENAME = goodrawdir + wave[i] + fn
    ;     ENDFOR

    ; ENDFOR

    ; ====================
    ; For correcting emcube satmaps.
    ; ====================

    satcube = MAKE_ARRAY(1001,1001,6)

    FOR i=0, N_ELEMENTS(emfiles)-1 DO BEGIN $       ; for each timestamp

        emfile = emfiles[i]                 
        r0 = STRPOS(emfile, '/', /reverse_search)
        fn = emfile.substring(r0 + 1)
        print, "TIMESTAMP: ", fn
        RESTORE, emfile

        FOR j=0, N_ELEMENTS(wave)-1 DO BEGIN        ; for each wavelength

            files = FILE_SEARCH(satdir + wave[j] + '/*.sav')      ; get the files for the wavelength
            file = files[i]                                       ; get proper timestamp
            r0 = STRPOS(files[i], '/', /reverse_search)                 
            print, "Getting ", wave[j], file.substring(r0)

            RESTORE, file         

            satcube[*,*,j] = satmap
            help, satmap, satcube

        ENDFOR
        print, ""
        satmap = satcube[x0:x1,y0:y1,*]
        help, satmap
        print, ""
        SAVE,  emcube, statuscube, satmap, logt, datacube, waveout, exptimeout, FILENAME = goodemdir + fn

    ENDFOR


END
    ; ; align satmap for each wavelength
    ; FOR j=0, (SIZE(datacube))[3]-1 DO BEGIN

    ;     ; but only if there's satmap to align!
    ;     where = WHERE(satmap[*,*,j], wherect)
    ;     IF wherect NE 0 THEN BEGIN

    ;         ; help, datacube
    ;         ; help, satmap
    ;         nx = N_ELEMENTS(datacube[*,0,0])
    ;         ny = N_ELEMENTS(datacube[0,*,0])

    ;         dc = FLTARR(nx,ny,2)
    ;         dc[*,*,1] = satmap[*,*,j]
    ;         dc[*,*,0] = datacube[*,*,j]
    ;         offset = tr_get_disp(dc, mad = 50, /debug)

    ;         newsatmap = image_translate(satmap[*,*,j], offset[*,1])
    ;         satcube[*,*,j] = newsatmap

    ;         print, "saturated"
    ;     ENDIF
    ; ENDFOR