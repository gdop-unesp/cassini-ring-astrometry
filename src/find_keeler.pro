pro find_keeler, et, ra, dec, pxscale, ang, distance, delta, pix
    keeler, et, 50, 0, 360, ring, lons
    cspice_reclat, ring.keelercen, radius, longs, lats
    ra = ra*cspice_rpd()
    dec = dec*cspice_rpd()
;    longs = longs/cspice_rpd()
;    lats = lats/cspice_rpd()
;    IF ra 	 GT !dpi THEN ra	-=2*!dpi
;    dra = longs - ra
;    dde = lats - dec
    dd = (acos(cos(dec)*cos(lats) + sin(dec)*sin(lats)*cos(ra-longs))/cspice_rpd())*60.0
;    dd = sqrt(dra*dra + dde*dde)
    k = where(dd EQ min(dd))
    distance = radius[k]
    ang = (lons[k])[0]/cspice_rpd()
    delt = (radius[k])[0]*tan(1024*pxscale)
    delta = atan(delt/136505.0)/cspice_rpd()
    pix = (radius[k])[0]*sin(pxscale)
    ;;print, longs[k], lats[k]
    
    return
   
end
