PRO brightpix

    dir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    files = FILE_SEARCH(dir + '*.sav')

    wave = ['94','131','171','193','211','335']

    FOR i=0, N_ELEMENTS(files)-1 DO BEGIN

        file = files[i]
        RESTORE, file
        print, file
        tot = TOTAL(datacube,3)

        im = IMAGE((total(emcube,3)^0.25), $
            ; POSITION = [0.55,0.1,0.95,0.95],$
            ; LAYOUT = [2,1,2], $
            /CURRENT)

        IF i NE 0 THEN BEGIN
            box.Delete
            ; c.Delete
        ENDIF

        ; c = CONTOUR(total(datacube, 3), $
        ;     C_LABEL_SHOW = 0,$
        ;     /OVERPLOT)

        ; c.getdata, brightpxls

        ; tot[WHERE(brightpxls EQ 0)] = 0

        brightest = MAX(tot, wheremax)
        help, tot
        print, "wheremax loc ", wheremax

        coords = ARRAY_INDICES(tot, wheremax)
        x = coords[0]
        y = coords[1]

        box = POLYGON([[x-10,y+10],[x+10,y+10],[x+10,y-10],[x-10,y-10]],$
            /DATA, FILL_BACKGROUND = 0, COLOR = 'red', TARGET = im) 


    ENDFOR


END