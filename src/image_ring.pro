pro imagem_view, file, dx=dx, dy=dy

    common SHARE, image

    IF not KEYWORD_SET(dx) THEN dx=0.0
    IF not KEYWORD_SET(dy) THEN dy=0.0

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
    rang = (acos(cos(dec)*cos(decb[0]) + sin(dec)*sin(decb[0])*cos(ra-rab[0]))/cspice_rpd())*60.0

    ;; get satelite position
    satpos, et, 'CASSINI', cod, pos
    ;; get saturn position
    satpos, et, 'CASSINI', '699', saturn

    ;; get stars positions
    catalogue, ra/cspice_rpd(), dec/cspice_rpd(), 3*rang, et, 'CASSINI', stars

    ;; get position and velocity of Cassini relative to Solar System Barycenter
    cspice_spkezr, 'CASSINI', et, 'J2000', 'NONE', '0', pos1, ltim

    ;; get keeler gap
    find_keeler, et, ra/cspice_rpd(), dec/cspice_rpd(), pxscale, ang, distance, deltax, interv ;; find the location of the keeler gap of the center of the image
    keeler, et, interv/5, ang-delta, ang+delta, ring, lons ;; get the keeler gap array
    refgap = where((lons - ang*cspice_rpd())^2 EQ min((lons-ang*cspice_rpd())^2)) ;; center of gap in the image, reference for calculation
    cspice_reclat, ring.keelercen[*,refgap], radiusref, raref, decref
    radec2pix, raref, decref, sampleref, lineref, CORRECT_STELAB=pos1[3:*]

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
    cspice_reclat, stars[0:2,*], radiusta, rasta, decsta
    radec2pix, rasta, decsta, samplesta, linesta, CORRECT_STELAB=pos1[3:*]

    volta = findgen(360)*cspice_rpd()
    ;; plot the image
    atv, image.image
        
    atvplot, samplekout+dx, linekout+dy, color='red', thick=3.0 ;; outer part of the gap
    index = WHERE((samplekout GE -200) and (samplekout LE 1250) and (linekout GE -200) and (linekout LE 1250), /NULL)
    atvplot, samplekin+dx, linekin+dy, color='red', thick=3.0 ;; inner part of the gap
    index = WHERE((samplekin GE -200) and (samplekin LE 1250) and (linekin GE -200) and (linekin LE 1250), /NULL)
    atvplot, sampleeout+dx, lineeout+dy, color='red', thick=3.0 ;; outer part of the gap
    index = WHERE((sampleeout GE -200) and (sampleeout LE 1250) and (lineeout GE -200) and (lineeout LE 1250), /NULL)
    atvplot, sampleein+dx, lineein+dy, color='red', thick=3.0 ;; inner part of the gap
    index = WHERE((sampleein GE -200) and (sampleein LE 1250) and (lineein GE -200) and (lineein LE 1250), /NULL)
    atvplot, samplearin+dx, linearin+dy, color='red', thick=3.0 ;; inner part of the gap
    index = WHERE((samplearin GE -200) and (samplearin LE 1250) and (linearin GE -200) and (linearin LE 1250), /NULL)
    FOR i=0, n_elements(samplesat) -1 do begin
        rr = atan(sat_rad[i]/radiusat[i])/image.pxscale
        atvplot, rr*sin(volta)+samplesat[i]+dx, rr*cos(volta)+linesat[i]+dy, color='red', thick=3.0
        atvxyouts, samplesat[i] + 5+dx, linesat[i] + 5+dy, sats[i], color='red', charsize=2
    ENDFOR
    foreach i, findgen(n_elements(samplesta)) do atvplot, 3.0*sin(volta)+samplesta[i]+dx, 3.0*cos(volta)+linesta[i]+dy, color='blue'

    cspice_kclear
end