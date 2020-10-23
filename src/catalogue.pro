pro catalogue, ra, dec, rang, et, obs, stars
    deg2mas = 1000.0d0*3600.0d0
    mas2deg = 1.0d0/deg2mas
    
    maglim = 13.0
    text = string(rang, maglim, ra, dec, format="vizquery -c.bm=%f -out=RA_ICRS -out=DE_ICRS -out=Plx -out=pmRA -out=pmDE -out=RV -out=Gmag -out=Epoch -source=I/345/gaia2 Gmag='<%4.1f' -c='%16.11f %16.11f' > cds.cat")
    print, 'Downloading catalogue stars'
    print, text
    spawn, text
    spawn, "egrep '(^0|^1|^2|^3|^4|^5|^6|^7|^8|^9|^ 0|^ 1|^ 2|^ 3|^ 4|^ 5|^ 6|^ 7|^ 8|^ 9|^  0|^  1|^  2|^  3|^  4|^  5|^  6|^  7|^  8|^  9)' cds.cat > gaia.cat"
    ;;read file
    openr, 1, 'gaia.cat'
    print, 'File "gaia.cat" generated'
    x =  [] & y = [] & z =  [] & mage = []
    i = 0
    cspice_spkpos, '0', et, 'J2000', 'NONE', obs, pos1, ltim
    print, 'Reading catalogue stars'
    while (~ EOF(1)) do begin
        ras = 0.0d0
        decs = 0.0d0
        plxs = 0.0001d0
        pmras = 0.0d0
        pmdes = 0.0d0
        rvs = 0.0d0
        epo = 0.0d0
        distance = 0.0d0
        readf, 1, $
        format = '(F15.11,1X,F15.11,1X,F10.4,1X,F9.3,1X,F9.3,1X,F7.2,1X,F7.4,1X,F6.1)', $
         ras, decs, plxs, pmras, pmdes, rvs, mags, epo
         
;        print, ras, decs, plxs, pmras, pmdes, rvs, mags, epo
        ;;Applying proper motion
        ras = ras*cspice_rpd() + (pmras/cos(decs*cspice_rpd()))*(et/cspice_jyear()-epo+2000.0d0)*mas2deg*cspice_rpd()
        decs = decs*cspice_rpd() + pmdes*(et/cspice_jyear()-epo+2000.0d0)*mas2deg*cspice_rpd()

        ;; Determining distance
        if (plxs lt 0.0001) then plxs = 0.0001d0
        distance = (1000.0d0/plxs)
        cspice_convrt, distance, "PARSECS", "KM", dista
;        print, distance, dista
        ;; Applying radial velocity
        dista = double(dista + rvs*(et - (epo-2000.0d0)*cspice_jyear()))
        
        cspice_latrec, dista, ras, decs, pos
        
;        print, ras/cspice_rpd(), decs/cspice_rpd(), dista, pos
         
        x = [x,pos[0] + pos1[0]]
        y = [y,pos[1] + pos1[1]]
        z = [z,pos[2] + pos1[2]]
        mage = [mage,mags]
        
        i = i+1
        
    endwhile
    print, 'Number of Gaia-DR2 stars read: ', i
    close, 1
    stars = [[x[0:i-1]],[y[0:i-1]],[z[0:i-1]],[mage[0:i-1]]]
    stars = transpose(stars)

end
