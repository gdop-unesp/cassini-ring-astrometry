pro compare_ring, lat1, lon1, lat2,lon2, dis, x, y
    x = []
    y = []
    ll = [[lon1], [lat1]]
    ll1 = reform(ll, 1024,1024,2)
    for i=0,1023 do begin
        ny = n_elements(lon2)
        x1 = transpose(ll1[i,*,0])#replicate(1,ny)
        y1 = lon2#replicate(1,1024)
        x2 = transpose(ll1[i,*,1])#replicate(1,ny)
        y2 = lat2#replicate(1,1024)
        k = where(abs(x1 - x2) lt step and abs(y1 - y2) lt step)
    endfor
    
end

