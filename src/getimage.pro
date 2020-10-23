pro getimage, file, image

    image = DICTIONARY()

    openr, arq, file,/get_lun
    a = assoc(arq, bytarr(4096))
    b = a[0]
    h = string(b)
    header = h.split('  ')
    remove, n_elements(header) - 1, header
    array = []
    array2 = []
    foreach element, header do begin
        if (strpos(element, '=') gt 0) then begin
            e = element.split('=')
            array = [array, e[0]]
            array2 = [array2, e[1]]
        endif else begin
            array2[-1] = array2[-1] + element
        endelse
    endforeach
    vars = hash(array,array2)
    image.header = vars

     ;; verifying camera information
    if (image.header['INSTRUMENT_ID'] eq "'ISSNA'") then begin
        image.code = -82360
        image.frame = "CASSINI_ISS_NAC"
        image.focal = 2002.703d0
        image.kmat = [[83.33333D, 0D, 0D], [0D, 83.3428D, 0D]]
        image.epsilon = [0D, 8.28D-6, 0D, 0D, 5.45D-6, -19.67D-6]
        image.center = [511.5d0, 511.5d0]
        image.bin = 1024/image.header['NL']
        imagem=fltarr(image.header['NL'],image.header['NL'])
        image.pxscale = (1.23d0/3600.0d0)*cspice_rpd()*image.bin
    endif else if (image.header['INSTRUMENT_ID'] eq "'ISSWA'") then begin
        image.code = -82361
        image.frame = "CASSINI_ISS_WAC"
        image.focal = 200.7761d0
        image.kmat = [[83.33333D, 0D, 0D], [0D, 83.34114D, 0D]]
        image.epsilon = [0D, 60.89D-6, 0D, 0D, 4.93D-6, -72.28D-6]
        image.center = [511.5d0, 511.5d0]
        image.bin = 1024/image.header['NL']
        imagem=fltarr(image.header['NL'],image.header['NL'])
        image.pxscale = (12.3d0/3600.0d0)*cspice_rpd()*image.bin
    endif

    point_lun,arq,4096+image.header['NBB']*4

    if image.header['FORMAT'] eq "'REAL'" then begin
        imagem=fltarr(fix(image.header['NL'])+image.header['NBB'],image.header['NL'])
        readu,arq,imagem
        image.image = imagem(image.header['NBB']:*,*)
    endif else if image.header['FORMAT'] eq "'BYTE'" then begin
        imagem=bytarr(fix(image.header['NL'])+ image.header['NBB'],image.header['NL'])
        readu,arq,imagem
        image.image = fix(imagem(image.header['NBB']:*,*))
    endif

;    point_lun,arq,4096
;    readu,arq,imagem
;    image.image = imagem

    ;imagem=bytarr(1024+24,1024)
    ;point_lun,arq, 4192
    ;readu,arq,imagem
    ;im=fix(imagem(24:*,*))

    close, arq

end
