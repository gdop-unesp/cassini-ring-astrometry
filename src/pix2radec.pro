PRO pix2radec_iter, sample, line, ra, dec, CORRECT_STELAB=vobs, VERBOSE=VERBOSE

    common SHARE, image

    pix2radec, sample, line, ra, dec

    nval = N_ELEMENTS(sample)
        
    grid_pitch = image.pxscale*image.bin
    grid_size = 50.0
    npts = grid_size*grid_size
    ra_gridVec = DBLARR(npts)
    dec_gridVec = DBLARR(npts)
    
    nite=3 
    
    
    FOR i=0, nval-1 DO BEGIN
        
        FOR k=0, nite-1 DO BEGIN
        
            IF KEYWORD_SET(VERBOSE) THEN PRINT, FORMAT='("Iteration #", I-, $)', k+1
            n_ite = 0
            REPEAT BEGIN
                n_ite++
            
                ra_vec  = ra[i]  + (INDGEN(grid_size)-grid_size/2) * grid_pitch
                dec_vec = dec[i] + (INDGEN(grid_size)-grid_size/2) * grid_pitch
                FOR j=0, grid_size-1 DO ra_gridVec[j*grid_size] = ra_vec
                FOR j=0, grid_size-1 DO dec_gridVec[j*grid_size] = REPLICATE(dec_vec[j], grid_size)
                        
                radec2pix, ra_gridVec, dec_gridVec, sample_gridVec, line_gridVec, CORRECT_STELAB=vobs
                residuals = (sample[i]-sample_gridVec)*(sample[i]-sample_gridVec) + (line[i]-line_gridVec)*(line[i]-line_gridVec)
                
                min_val = MIN(residuals, index)
                min_pos = ARRAY_INDICES(INTARR(grid_size, grid_size, /NOZERO), index)
                ra[i]  = ra_vec[min_pos[0]]
                dec[i] = dec_vec[min_pos[1]]
            
                IF KEYWORD_SET(VERBOSE) THEN BEGIN
                    format = '(" Min residual (pixels): ", G-12.6, " at (", I3, ",", I3, ")")'
                    PRINT, FORMAT=format, min_val, min_pos[0], min_pos[1]
                ENDIF
            ENDREP UNTIL min_pos[0] GT 0 || min_pos[0] LT grid_size-1 || $
                         min_pos[1] GT 0 || min_pos[1] LT grid_size-1 || $
                         n_ite GE 10
            IF n_ite EQ 10 THEN break
        
            grid_pitch *= 0.1D
        ENDFOR
        
    ENDFOR
    
    RETURN
END


pro pix2radec, sample, line, ra, dec

    common SHARE, image

    cspice_invert, image.matinv, mat

	n = N_ELEMENTS(sample)

	K = image.kmat * (DOUBLE(image.focal)/image.bin)
    kt = transpose(image.kmat)
    invK = TOTAL(image.kmat[2,*]) EQ 0 ? LA_INVERT(image.kmat[0:1,*]) : LA_INVERT(Kt##K)##Kt
    
    xy = invK ## [[-sample+image.center[0]], [-line+image.center[1]]]
    
    rhoCAM = MAKE_ARRAY(n, 3, VALUE=1.0D)
	rhoCAM[*,0:1] = xy[*,0:1]/image.focal
	
	invNorm_rho_cam = DBLARR(n,3)
	FOR i=0, n-1 DO invNorm_rho_cam[i, *] = 1.0D/NORM(rhoCAM[i,*])
	
	rhoJ2000 = mat##(rhoCAM*invNorm_rho_cam)
	
	cspice_recrad, TRANSPOSE(rhoJ2000), range, ra, dec
	
    return
end

