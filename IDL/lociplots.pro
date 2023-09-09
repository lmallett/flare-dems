PRO lociplots, coords = coords

    IF NOT KEYWORD_SET(coords) THEN BEGIN
        x = 125
        y = 275
        coords = [x,y]
    ENDIF ELSE BEGIN
        x = coords[0]
        y = coords[1]
    ENDELSE

    print, coords

    dir = '/home/reu23-student/Desktop/Research/2011_08_09/emcubes/all/'
    files = FILE_SEARCH(dir + '*.sav')
    
    aia_tresp = make_aia_response('20110809',/feldman)
    help, aia_tresp

    wave = ['94','131','171','193','211','335']
    colors = ['r','g','b', 'c', 'm', 'y']

    w = WINDOW(DIMENSIONS = [1280,720])

    FOR i=0, N_ELEMENTS(files)-1 DO BEGIN

        file = files[i]
        RESTORE, file
        print, i
        print, file
        ; CATCH, error_status
        ; ; write error handler
        ; IF emcube[coords[0]] 
        ; logt2 = cgNumber_Formatter(logt, decimals = 2)

        ; IF emcube[x,y,0] EQ 0 THEN CONTINUE

        p0 = PLOT(logt, emcube[x,y,*], $
            XMINOR = 0, $
            XTICKNAME = logt, $
            XTICKVALUES = logt, $
            XTICKFORMAT = '(F10.1)', $
            XRANGE = [min(logt), max(logt)],$
            XTITLE = '$log(T)$', $
            /YLOG,$
            YRANGE = [1,1000000],$
            YSTYLE = 1,$
            YTITLE = 'Emission Measure $[10^{26} cm^{-5}]$', $
            /HISTOGRAM, $
            POSITION = [0.1,0.1,0.55,0.95],$
            ; DIMENSIONS = [1280,720],$
            LAYOUT = [2,1,1], $
            ; MARGIN = [0.2,0.1,0,0.05], $ ; left bottom right top
            NAME = 'EM', /CURRENT)

        wave = ['94','131','171','193','211','335']
        colors = ['r','g','b', 'c', 'm', 'y']

        FOR j=0, N_ELEMENTS(wave)-1 DO BEGIN

            response = 'aia_tresp.resp_' + wave[j]
            print, response
            datarray = MAKE_ARRAY(20, /FLOAT, value = datacube[x,y,j])
            
            void = EXECUTE('datarray2 = datarray / float(' + response + ')')
            datarray2 = datarray2/(10.^(26))
            
            void2  = EXECUTE('p' + strtrim(j+1,2) + "= PLOT(logt, datarray2," +$
                "COLOR = colors[j], LAYOUT =[2,1,1],"+$
                "/OVERPLOT, NAME ='"+wave[j]+"')")

        ENDFOR
        timestamp = file.substring(-27,-5)

        leg = LEGEND(TRANSPARENCY = 0,$
            ; POSITION = [0.52,0.9],$
            ; POSITION = [0.23,0.4],$
            POSITION = [260,650],$
            SAMPLE_WIDTH = 0.1,$
            ; ORIENTATION = 1,$
            /DEVICE)
        ; text01 = TEXT(0.12,0.87, timestamp, /FILL_BACKGROUND)

        im = IMAGE(bytscl(total(emcube,3)^0.25), $
            POSITION = [0.55,0.1,0.95,0.95],$
            ; MARGIN = [0, 0.1, 0, 0.05], $
            LAYOUT = [2,1,2], /CURRENT)

        box = POLYGON([[x-10,y+10],[x+10,y+10],[x+10,y-10],[x-10,y-10]],$
            /DATA, FILL_BACKGROUND = 0, COLOR = 'red', TARGET = im) 

        ; text01 = TEXT(300,700, timestamp, /FILL_BACKGROUND, /DEVICE)
        ; text02print = TEXT(0.04,0.025, text02, COLOR = [255,255,255], TARGET=im, /RELATIVE)
        ; text01 = TEXT(0.04, 0.05, timestamp, COLOR = [255,255,255], TARGET = im, /RELATIVE)

        text01 = "[" + STRTRIM(coords[0],2) + "," + STRTRIM(coords[1],2) + "]"
        text01 = timestamp + "$\n$" + text01
        text01print = TEXT(0.04, 0.03, text01, COLOR = [255,255,255], TARGET = im, /RELATIVE)

        p0.Save, "/home/reu23-student/Desktop/Research/2011_08_09/lociplots140x240/" + timestamp + ".png"
        p0.Erase

        ; PLOT, logt, emcube[x,y,*], $
        ;     xtickname = logt, $
        ;     xtickv = logt, $
        ;     xrange = [min(logt), max(logt)],$
        ;     xtickformat = '(F10.1)', $
        ;     psym = 10
        ; loadct, 8
        ; WINDOW, 30, XSIZE = (SIZE(emcube))[1], YSIZE = (SIZE(emcube))[2]
        ; TV, BYTSCL(TOTAL(emcube,3)^0.25)
        ; POLYFILL, [x-20, x-20, x+20, x+100], [y-20, y+20, y-20, y+20], /LINE_FILL, COLOR =0, /DATA
        ; POLYFILL, [0.1, 0.2, 0.1, 0.2], [0.1, 0.2, 0.2, 0.1], /LINE_FILL, COLOR =0
        ; DRAW_ROI, [0.1, 0.2, 0.1, 0.2], [0.1, 0.2, 0.2, 0.1], /LINE_FILL, COLOR =0, /DATA
    ENDFOR    
END