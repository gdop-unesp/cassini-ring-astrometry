pro keeler, et, step, angmin, angmax, ring, lons

    ring = DICTIONARY()

    au2km = 1.495978707d8
    day2s = 86400.0d0
    s2day = 1.0/day2s
    ltvel = 1.7314463267424034d2*(au2km/day2s)
    rkeeler = 136502.0 ;; km
    wkeeler = 37.0 ;; km
    rencke = 133583 ;; km
    wencke = 321 ;; km
    aring = 136770 ;; km
    lgap = 2.0*!DPI*rkeeler ;; lenght of the gap
    int = 2.0*!DPI*step/lgap ;; step in angle
    if angmax LT angmin then angmin = angmin - 360.
    angmin = angmin*cspice_rpd()
    angmax = angmax*cspice_rpd()
    n = (angmax-angmin)/int ;; number of intervals within the limits given.

    T = et/(day2s*36525.0)
    d = et*s2day
    alfa0 = 40.589 - 0.036*T
    delta0 = 83.537 - 0.004*T
    W = 38.9 + 810.7939024*d

    polea = alfa0*cspice_rpd()
    poled = delta0*cspice_rpd()
    cspice_rotate, (poled - 90.0*cspice_rpd())[0], 1, mx
    cspice_rotate, (- polea - 90.0*cspice_rpd())[0], 3, mz
    cspice_pxform, "IAU_SATURN", "J2000", et, matr
    
    lons = findgen(n)*(angmax-angmin)/n + angmin
    lats = replicate(0.0,n)
    cspice_latrec, replicate(rkeeler,n), lons, lats, ckee
    cspice_latrec, replicate(rkeeler - wkeeler/2.0,n), lons, lats, ikee
    cspice_latrec, replicate(rkeeler + wkeeler/2.0,n), lons, lats, okee
    cspice_latrec, replicate(rencke,n), lons, lats, cenc
    cspice_latrec, replicate(rencke - wencke/2.0,n), lons, lats, ienc
    cspice_latrec, replicate(rencke + wencke/2.0,n), lons, lats, oenc
    cspice_latrec, replicate(aring,n), lons, lats, oarin
;;    cent = transpose(mx##transpose(cent))
;;    cent = transpose(mz##transpose(cent))
    ckee = transpose(matr##transpose(ckee))
;;    bai = transpose(mx##transpose(bai))
;;    bai = transpose(mz##transpose(bai))
    ikee = transpose(matr##transpose(ikee))
;;    cim = transpose(mx##transpose(cim))
;;    cim = transpose(mz##transpose(cim))
    okee = transpose(matr##transpose(okee))
    cenc = transpose(matr##transpose(cenc))
    ienc = transpose(matr##transpose(ienc))
    oenc = transpose(matr##transpose(oenc))
    oarin = transpose(matr##transpose(oarin))
    
    eta = replicate(et, n)

    cspice_spkpos, '0', eta, 'J2000', 'NONE', 'Cassini', pos1, ltim

    ltime = replicate(0.0d0, n)
    
    tt = 0
    
    ;; loop
    repeat begin
    
        etn = eta - ltime
        
        cspice_spkpos, 'Saturn', etn, 'J2000', 'NONE', '0', pos2, ltim
        
        ring.keelercen = pos1 + pos2 + ckee
        ring.keelerin = pos1 + pos2 + ikee
        ring.keelerout = pos1 + pos2 + okee
        ring.enckecen = pos1 + pos2 + cenc
        ring.enckein = pos1 + pos2 + ienc
        ring.enckeout = pos1 + pos2 + oenc
        ring.aring = pos1 + pos2 + oarin
        
        ltime = sqrt(ring.keelercen[0,*]^2+ring.keelercen[1,*]^2+ring.keelercen[2,*]^2)/ltvel
    
    endrep until array_equal(abs(eta - etn - ltime) LT 0.001, 1)  EQ 1
    ;; fim loop
    
    return

end