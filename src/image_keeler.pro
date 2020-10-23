pro roda, file, dx=dx, dy=dy, rangex=rangex, rangey=rangey, incrx=incrx, incry=incry, output=matrix_sum, high=high, low=low, sigma=sigma, ARING=aring, OKEELER=okeeler, IKEELER=ikeeler, OENCKE=oencke, IENCKE=iencke, CATAG=catag

    common SHARE, image
    
    IF not KEYWORD_SET(dx) THEN dx=0.0
    IF not KEYWORD_SET(dy) THEN dy=0.0
    IF not KEYWORD_SET(rangex) AND not KEYWORD_SET(rangey) THEN BEGIN
        rangex=250
        rangey=250
    ENDIF
    IF not KEYWORD_SET(incrx) AND not KEYWORD_SET(incry) THEN BEGIN
        incrx=0.5
        incry=0.5
    ENDIF
    IF not KEYWORD_SET(high) THEN high=0.985
    IF not KEYWORD_SET(low) THEN low=0.977
    IF not KEYWORD_SET(sigma) THEN sigma=0.05

    delta = 90.0
    cod =     ['618', '635',     '615',   '616',        '617',     '611',        '610',   '653',     '601',   '632',     '649',   '633',     '602',       '603',    '613',     '614',     '634',        '604',   '612',    '605',  '606',   '607',      '608',     '609']
    sats =    ['Pan', 'Daphnis', 'Atlas', 'Prometheus', 'Pandora', 'Epimetheus', 'Janus', 'Aegaeon', 'Mimas', 'Methone', 'Anthe', 'Pallene', 'Enceladus', 'Tethys', 'Telesto', 'Calypso', 'Polydeuces', 'Dione', 'Helene', 'Rhea', 'Titan', 'Hyperion', 'Iapetus', 'Phoebe']
    sat_rad = [14.1,   3.8,       15.1,    43.1,         40.7,      58.1,         89.5,    0.3,       198.2,   1.6,       0.9,     2.5,       252.1,       533.0,    12.4,      10.7,      1.3,          561.7,   17.6,     764.3,  2574.7,  135.0,      735.6,     106.5]

    ;; open image and get header
    getimage, file, image

    ;; open cspice and get dates
    cspice_furnsh, '/home/altair/idl/caviar_r1.0/kernels.ker' ;; open kernels
    utc1 = image.header['START_TIME'].replace('Z', '') ;; get image instant
    utc1 = utc1.replace("'", '')
    utc2 = image.header['STOP_TIME'].replace('Z', '') ;; get image instant
    utc2 = utc2.replace("'", '')
    cspice_str2et, utc1, et1 ;; convert utc to et
    cspice_str2et, utc2, et2 ;; convert utc to et
    et = (et1 + et2)/2.0
    image.et = et
    cspice_et2utc, et, 'ISOC', 3, utc
    print, 'Central instant: ', utc

    ;; get FOV
    cspice_getfov, image.code, 4, shape, frame, bsight, bounds ;; get fov of cassini
    pxscale = bounds[0]/512 ;; pixel scale
    cspice_pxform, "J2000", image.frame, et, matinv ;; matrix icrs -> cassini
    cspice_m2eul, matinv, 3, 1, 3, twist, dec, ra

    ;; correct twist angle from tajeddine (2013) = -9.6d-2 degree
    dtwist = 0.0
    image.twist = twist + dtwist*cspice_rpd()
    cspice_eul2m, image.twist, dec, ra, 3, 1, 3, matinv
    cspice_invert, matinv, mat
    matini = matinv
    dec   = !dpi/2 - dec
    ra    = ra - !dpi/2
    IF ra    LT 0.0d0 THEN ra   +=2*!dpi
    image.ra = ra
    image.dec = dec
    image.twist = twist
    image.matinv = matinv
    image.mat = mat
    print, 'Initial center RA=', image.ra/cspice_rpd(), ' deg, DEC=', image.dec/cspice_rpd(), ' deg'
    print, 'Initial center RA=', image.ra, ' rad, DEC=', image.dec, ' rad'
    print, 'Twist=', twist/cspice_rpd()
    coordbou = transpose(image.mat##transpose(bounds))
    cspice_reclat, coordbou, radb, rab, decb ;; get the bounds alfa and delta
    rang = (acos(cos(dec)*cos(decb[0]) + sin(dec)*sin(decb[0])*cos(ra-rab[0]))/cspice_rpd())*60.0 ;; angle from the corner to the center.

    ;; get satelite position
    satpos, et, 'CASSINI', cod, pos
    ;; get saturn position
    satpos, et, 'CASSINI', '699', saturn

    ;; get stars positions
    IF KEYWORD_SET(CATAG) THEN BEGIN
        catalogue, ra/cspice_rpd(), dec/cspice_rpd(), 3*rang, et, 'CASSINI', stars
    ENDIF

    ;; get position and velocity of Cassini relative to Solar System Barycenter
    cspice_spkezr, 'CASSINI', et, 'J2000', 'NONE', '0', pos1, ltim

    ;; get keeler gap
    find_keeler, et, ra/cspice_rpd(), dec/cspice_rpd(), pxscale, ang, distance, deltax, interv ;; find the location of the keeler gap of the center of the image
    keeler, et, interv/5, ang-delta, ang+delta, ring, lons ;; get the keeler gap array
    refgap = where((lons - ang*cspice_rpd())^2 EQ min((lons-ang*cspice_rpd())^2)) ;; center of gap in the image, reference for calculation
    cspice_reclat, ring.keelercen[*,refgap], radiusref, raref, decref
    radec2pix, raref, decref, sampleref, lineref, CORRECT_STELAB=pos1[3:*]

    print, radiusref*sin(image.pxscale), ' km'

    deltax = findgen(rangex/incrx, increment=incrx) - rangex/2 + dx
    deltay = findgen(rangey/incry, increment=incry) - rangey/2 + dy
    
    print, 'x = [', min(deltax), ',', max(deltax), '], y = [', min(deltay), ',', max(deltay), ']'
;    deltax = [-23]
;    deltay = [-1]
    matrix_sum = MAKE_ARRAY(n_elements(deltax), n_elements(deltay), /DOUBLE, VALUE = 0)

    radec2pix, ra, dec, centerx, centery
    
    ;; get image with the canny filter
    filtered = CANNY(image.image, high=high, low=low, sigma=sigma)
    ;atv, filtered
    ;atvplot, samplecen, linecen, color='red' ;; center of the gap
    ;atvplot, samplecim, linecim, color='blue' ;; outer part of the gap
    ;atvplot, samplebai, linebai, color='yellow' ;; inner part of the gap

     ;; calculates position on the image from the center of keeler gap
    cspice_reclat, ring.keelercen, radiuscen, rakcen, deckcen
    radec2pix, rakcen, deckcen, samplekcen, linekcen, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the inner part of the keeler gap
    cspice_reclat, ring.keelerin, radiusbai, rakin, deckin
    radec2pix, rakin, deckin, samplekin, linekin, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the outer part of the keeler gap
    cspice_reclat, ring.keelerout, radiuscim, rakout, deckout
    radec2pix, rakout, deckout, samplekout, linekout, CORRECT_STELAB=pos1[3:*]
    ;; calculates position on the image from the center of keeler gap
    cspice_reclat, ring.enckecen, radiuscen, raecen, dececen
    radec2pix, raecen, dececen, sampleecen, lineecen, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the inner part of the keeler gap
    cspice_reclat, ring.enckein, radiusbai, raein, decein
    radec2pix, raein, decein, sampleein, lineein, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the outer part of the keeler gap
    cspice_reclat, ring.enckeout, radiuscim, raeout, deceout
    radec2pix, raeout, deceout, sampleeout, lineeout, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the outer part of the A ring
    cspice_reclat, ring.aring, radiusarin, raarin, decarin
    radec2pix, raarin, decarin, samplearin, linearin, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the satellite
    cspice_reclat, pos[0:2,*], radiusat, rasat, decsat
    radec2pix, rasat, decsat, samplesat, linesat, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the catalogue stars
    IF KEYWORD_SET(CATAG) THEN BEGIN
        cspice_reclat, stars[0:2,*], radiusta, rasta, decsta
        radec2pix, rasta, decsta, samplesta, linesta, CORRECT_STELAB=pos1[3:*]
    ENDIF
            
    cgProgressBar = Obj_New("cgProgressBar")
    cgProgressBar -> Start
    nn = n_elements(deltax)*n_elements(deltay)
    k = 0
    FOR i = 0, n_elements(deltax)-1 do begin
        dx = deltax[i]
        FOR j = 0, n_elements(deltay)-1 DO begin

            dy = deltay[j]

            ; create a image where the calculated rings are 1 and the rest are 0
            teste = MAKE_ARRAY(1024, 1024, /INTEGER, VALUE = 0)
            IF KEYWORD_SET(ARING) THEN BEGIN
                sa = round(samplearin + dx)
                la = round(linearin + dy)
                index = WHERE((sa GE 0) and (sa LE 1023) and (la GE 0) and (la LE 1023), /NULL)
                teste[sa[index],la[index]] = 1
            ENDIF
            IF KEYWORD_SET(IKEELER) THEN BEGIN
                sb = round(samplekin + dx)
                lb = round(linekin + dy)
                index = WHERE((sb GE 0) and (sb LE 1023) and (lb GE 0) and (lb LE 1023), /NULL)
                teste[sb[index],lb[index]] = 1
            ENDIF
            IF KEYWORD_SET(OKEELER) THEN BEGIN
                sc = round(samplekout + dx)
                lc = round(linekout + dy)
                index = WHERE((sc GE 0) and (sc LE 1023) and (lc GE 0) and (lc LE 1023), /NULL)
                teste[sc[index],lc[index]] = 1
            ENDIF
            IF KEYWORD_SET(IENCKE) THEN BEGIN
                sd = round(sampleein + dx)
                ld = round(lineein + dy)
                index = WHERE((sd GE 0) and (sd LE 1023) and (ld GE 0) and (ld LE 1023), /NULL)
                teste[sd[index],ld[index]] = 1
            ENDIF
            IF KEYWORD_SET(OENCKE) THEN BEGIN
                se = round(sampleeout + dx)
                lee = round(lineeout + dy)
                index = WHERE((se GE 0) and (se LE 1023) and (lee GE 0) and (lee LE 1023), /NULL)
                teste[se[index],lee[index]] = 1
            ENDIF

            new_frame = teste*filtered
            
            matrix_sum[i,j] = double(total(new_frame))

            k = k + 1
            
            cgProgressBar -> Update, (double(k)/double(nn))*100
;            print, i, j, dx, dy, total(new_frame), total(teste)

        ENDFOR
    ENDFOR
    cgProgressBar -> Destroy

    mx = MAX(matrix_sum, location)
    loc = ARRAY_INDICES(matrix_sum, location)
    print, loc

    ind = where(matrix_sum/mx gt 0.75)
    s = SIZE(matrix_sum)
    ncol = s(1)
    col = ind MOD ncol
    row = ind / ncol

    xcen = TOTAL((col*matrix_sum[ind]))/TOTAL(matrix_sum[ind])
    ycen = TOTAL((row*matrix_sum[ind]))/TOTAL(matrix_sum[ind])

    print, xcen, ycen

    dx = deltax[xcen]
    dy = deltay[ycen]
    print, "Correction in pixels ", dx, dy

    ;pars = []

    ;mask = MAKE_ARRAY(n_elements(deltax), n_elements(deltay), VALUE = 0)
    ;mask[ind] = 1.0

    ;GAUSS2DFIT( matrix_sum, pars, MASK=mask , /TILT )
    ;print, pars
    ;print, xcen, ycen


    OPENW, 4, 'variaveis.py'
    PRINTF, 4, 'import numpy as np'

    image.matinv = matini

    volta = findgen(360)*cspice_rpd()
    ;; plot the image
    atv, image.image ;; image
;    atv, filtered
    ;atvplot, samplecen, linecen, color='red' ;; center of the gap
    IF KEYWORD_SET(OKEELER) THEN begin 
        atvplot, samplekout, linekout, color='red', thick=3.0 ;; outer part of the gap
        PRINTF, 4, 'inkout = np.array(['
        index = WHERE((samplekout GE -200) and (samplekout LE 1250) and (linekout GE -200) and (linekout LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', samplekout[i], ',', linekout[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(IKEELER) THEN begin 
        atvplot, samplekin, linekin, color='red', thick=3.0 ;; inner part of the gap
        PRINTF, 4, 'inkin = np.array(['
        index = WHERE((samplekin GE -200) and (samplekin LE 1250) and (linekin GE -200) and (linekin LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', samplekin[i], ',', linekin[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(OENCKE) THEN begin 
        atvplot, sampleeout, lineeout, color='red', thick=3.0 ;; outer part of the gap
        PRINTF, 4, 'ineout = np.array(['
        index = WHERE((sampleeout GE -200) and (sampleeout LE 1250) and (lineeout GE -200) and (lineeout LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', sampleeout[i], ',', lineeout[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(IENCKE) THEN begin 
        atvplot, sampleein, lineein, color='red', thick=3.0 ;; inner part of the gap
        PRINTF, 4, 'inein = np.array(['
        index = WHERE((sampleein GE -200) and (sampleein LE 1250) and (lineein GE -200) and (lineein LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', sampleein[i], ',', lineein[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(ARING) THEN begin 
        atvplot, samplearin, linearin, color='red', thick=3.0 ;; inner part of the gap
        PRINTF, 4, 'inarin = np.array(['
        index = WHERE((samplearin GE -200) and (samplearin LE 1250) and (linearin GE -200) and (linearin LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', samplearin[i], ',', linearin[i], '],'
        PRINTF, 4, '])'
    endif
    PRINTF, 4, 'insats = np.array(['
    FOR i=0, n_elements(samplesat) -1 do begin
        rr = atan(sat_rad[i]/radiusat[i])/image.pxscale
        atvplot, rr*sin(volta)+samplesat[i], rr*cos(volta)+linesat[i], color='red', thick=3.0
        PRINTF, 4, '[', samplesat[i], ',', linesat[i], '],'
;        if (samplesat[i] GE 0) and (samplesat[i] LE 1023) and (linesat[i] GE 0) and (linesat[i] LE 1023) then begin
        ;atvxyouts, samplesat[i] + 5, linesat[i] + 5, sats[i], color='yellow', charsize=2
;        endif
    ENDFOR
    PRINTF, 4, '])'
    ;atvplot, samplecen, linecen, color='yellow' ;; center of the gap

    pix2radec_iter, image.center[0] - dx, image.center[1] - dy, new_ra, new_dec
    cspice_eul2m, image.twist, 0.5*!dpi-new_dec, 0.5*!dpi+new_ra, 3, 1, 3, matinv
    cspice_invert, matinv, mat
    IF new_ra LT 0.0d0 THEN new_ra  +=2*!dpi
    ;IF twist LT 0.0d0 THEN twist+=2*!dpi
    print, 'center corrected RA=', (new_ra/cspice_rpd()), ' deg, DEC=', new_dec/cspice_rpd(), ' deg'
    print, 'center corrected RA=', new_ra, ' rad, DEC=', new_dec, ' rad'
    print, 'delta position dRA =', ((new_ra-ra)/cspice_rpd())*3600.0d0, ' arcsec DEC=', ((new_dec-dec)/cspice_rpd())*3600.0d0, ' arcsec'
    ;image.ra = new_ra
    ;image.dec = new_dec
    image.matinv = matinv

    ;;;;;;; Calculating with new center

     ;; calculates position on the image from the center of keeler gap
    cspice_reclat, ring.keelercen, radiuscen, rakcen, deckcen
    radec2pix, rakcen, deckcen, samplekcen, linekcen, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the inner part of the keeler gap
    cspice_reclat, ring.keelerin, radiusbai, rakin, deckin
    radec2pix, rakin, deckin, samplekin, linekin, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the outer part of the keeler gap
    cspice_reclat, ring.keelerout, radiuscim, rakout, deckout
    radec2pix, rakout, deckout, samplekout, linekout, CORRECT_STELAB=pos1[3:*]
    ;; calculates position on the image from the center of keeler gap
    cspice_reclat, ring.enckecen, radiuscen, raecen, dececen
    radec2pix, raecen, dececen, sampleecen, lineecen, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the inner part of the keeler gap
    cspice_reclat, ring.enckein, radiusbai, raein, decein
    radec2pix, raein, decein, sampleein, lineein, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the outer part of the keeler gap
    cspice_reclat, ring.enckeout, radiuscim, raeout, deceout
    radec2pix, raeout, deceout, sampleeout, lineeout, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the outer part of the A ring
    cspice_reclat, ring.aring, radiusarin, raarin, decarin
    radec2pix, raarin, decarin, samplearin, linearin, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the satellite
    cspice_reclat, pos[0:2,*], radiusat, rasat, decsat
    radec2pix, rasat, decsat, samplesat, linesat, CORRECT_STELAB=pos1[3:*]
    ;; does the same for the catalogue stars
    IF KEYWORD_SET(CATAG) THEN BEGIN
        cspice_reclat, stars[0:2,*], radiusta, rasta, decsta
        radec2pix, rasta, decsta, samplesta, linesta, CORRECT_STELAB=pos1[3:*]
    ENDIF

    IF KEYWORD_SET(OKEELER) THEN begin 
        atvplot, samplekout, linekout, color='yellow', thick=3.0 ;; outer part of the gap
        PRINTF, 4, 'fikout = np.array(['
        index = WHERE((samplekout GE -200) and (samplekout LE 1250) and (linekout GE -200) and (linekout LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', samplekout[i], ',', linekout[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(IKEELER) THEN begin
        atvplot, samplekin, linekin, color='yellow', thick=3.0 ;; inner part of the gap
        PRINTF, 4, 'fikin = np.array(['
        index = WHERE((samplekin GE -200) and (samplekin LE 1250) and (linekin GE -200) and (linekin LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', samplekin[i], ',', linekin[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(OENCKE) THEN begin
        atvplot, sampleeout, lineeout, color='yellow', thick=3.0 ;; outer part of the gap
        PRINTF, 4, 'fieout = np.array(['
        index = WHERE((sampleeout GE -200) and (sampleeout LE 1250) and (lineeout GE -200) and (lineeout LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', sampleeout[i], ',', lineeout[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(IENCKE) THEN begin
        atvplot, sampleein, lineein, color='yellow', thick=3.0 ;; inner part of the gap
        PRINTF, 4, 'fiein = np.array(['
        index = WHERE((sampleein GE -200) and (sampleein LE 1250) and (lineein GE -200) and (lineein LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', sampleein[i], ',', lineein[i], '],'
        PRINTF, 4, '])'
    endif
    IF KEYWORD_SET(ARING) THEN begin
        atvplot, samplearin, linearin, color='yellow', thick=3.0 ;; inner part of the gap
        PRINTF, 4, 'fiarin = np.array(['
        index = WHERE((samplearin GE -200) and (samplearin LE 1250) and (linearin GE -200) and (linearin LE 1250), /NULL)
        foreach i, index do PRINTF, 4, '[', samplearin[i], ',', linearin[i], '],'
        PRINTF, 4, '])'
    endif
    PRINTF, 4, 'fisats = np.array(['
    FOR i=0, n_elements(samplesat) -1 do begin
        rr = atan(sat_rad[i]/radiusat[i])/image.pxscale
        atvplot, rr*sin(volta)+samplesat[i], rr*cos(volta)+linesat[i], color='yellow', thick=3.0
        PRINTF, 4, '[', samplesat[i], ',', linesat[i], '],'
;        if (samplesat[i] GE 0) and (samplesat[i] LE 1023) and (linesat[i] GE 0) and (linesat[i] LE 1023) then begin
        atvxyouts, samplesat[i] + 5, linesat[i] + 5, sats[i], color='yellow', charsize=2
;        endif
    ENDFOR
    PRINTF, 4, '])'
    IF KEYWORD_SET(CATAG) THEN BEGIN
       foreach i, findgen(n_elements(samplesta)) do atvplot, 3.0*sin(volta)+samplesta[i], 3.0*cos(volta)+linesta[i], color='blue'
    ENDIF
    
    FREE_LUN, 4

    cspice_kclear

    writefits, 'matrix.fits', matrix_sum
    writefits, 'imagem.fits', image.image
    writefits, 'canny.fits', filtered
    
    return
    

end


pro check_canny, file, filtered, high=high, low=low, sigma=sigma
    IF not KEYWORD_SET(high) THEN high=0.985
    IF not KEYWORD_SET(low) THEN low=0.977
    IF not KEYWORD_SET(sigma) THEN sigma=0.05
    getimage, file, image
    filtered = CANNY(image.image, high=high, low=low, sigma=sigma)
end