pro satpos, et, obs, cod, pos
    au2km = 1.495978707d8
    day2s = 86400.0d0
    s2day = 1.0d0/day2s
    ltvel = 1.7314463267424034d2*(au2km/day2s)
    
    eta = replicate(et, n_elements(cod))

    cspice_spkezr, '0', et, 'J2000', 'NONE', obs, pos1, ltim

    ltime = replicate(0.0d0, n_elements(cod))
    
    tt = 0
    
    pos = double(fltarr(6,n_elements(cod)))

    print, 'Obtaining satellite positions'
    
    ;; loop
    repeat begin
    
        etn = eta - ltime
        
        FOREACH i, findgen(n_elements(cod)) DO BEGIN
            cspice_spkezr, cod[i], etn[i], 'J2000', 'NONE', '0', pos2, ltim
            pos[0,i] = pos1[0] + pos2[0]
            pos[1,i] = pos1[1] + pos2[1]
            pos[2,i] = pos1[2] + pos2[2]
            pos[3,i] = pos1[3] + pos2[3]
            pos[4,i] = pos1[4] + pos2[4]
            pos[5,i] = pos1[5] + pos2[5]
        ENDFOREACH
        
        ltime = sqrt(pos[0,*]^2+pos[1,*]^2+pos[2,*]^2)/ltvel
    
    endrep until array_equal(abs(eta - etn - ltime) LT 0.001, 1)  EQ 1
    ;; fim loop

end
