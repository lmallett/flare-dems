;+
; Gets datacubes from emcubes
;
; WRITTEN BY: Lucien Mallett
;-

PRO export_datacubes

    dir     = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    outdir  = '/home/reu23-student/Desktop/Research/2011_08_09/datacubes/'

    files = FILE_SEARCH(dir + '*.sav')
    FOREACH file, files DO BEGIN
        RESTORE, file
        file = file.substring(59)
        print, file
        SAVE, datacube, FILENAME = outdir + file
    ENDFOREACH

END