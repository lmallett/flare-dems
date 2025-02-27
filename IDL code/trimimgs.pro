;+
; Helper function. Crops a dataset. All parts of the dataset
; (datacube, emcube, statuscube, satmap) are assumed to be the
; same size (though this can be edited).
;
; WRITTEN BY: Lucien Mallett
;
; LAST MODIFIED: 2023-07-20
;-

PRO trimimgs

    ; xregion = [350,650]
    ; yregion = [275,775]
    ; dir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all-uncropped/'
    ; outdir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'

    xregion = [250,750]
    yregion = [350,750]    
    dir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all-uncropped/'
    outdir = '/home/reu23-student/Desktop/Research/2014_09_10/emcubes/all/'

    ;;;;;;
    xstart  = xregion[0]
    xend    = xregion[1]-1
    ystart  = yregion[0]
    yend    = yregion[1]-1

    files = FILE_SEARCH(dir + '*.sav')

    FOREACH file, files DO BEGIN

        RESTORE, file
        datacube = TEMPORARY(datacube[xstart:xend, ystart:yend, *])
        emcube = TEMPORARY(emcube[xstart:xend, ystart:yend, *])
        statuscube = TEMPORARY(statuscube[xstart:xend, ystart:yend, *])
        satmap = TEMPORARY(satmap[xstart:xend, ystart:yend, *])

        fn = file.substring(-38)
        print, outdir + fn
        SAVE, emcube, statuscube, satmap, logt, datacube, waveout, exptimeout, FILENAME = outdir + fn

    ENDFOREACH

END