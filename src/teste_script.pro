.r remove ;; function that remove elements from array
.r getimage ;; get image function
.r keeler ;; get keeler gap function
.r find_keeler ;; find keeler gap function
.r satpos ;; get position of a satellite
.r catalogue ;; get catalogue stars
.r radec2pix ;; ra, dec -> sample, line
.r pix2radec ;; sample, line -> ra, dec
.r atv ;; run atv plot aplication
.r image_keeler ;; code to fit a image
.r cgprogressbar__define ;; show progressbar
.r writefits ;; write fits file
.r image_ring

file = 'imagens/N1502843284_1_CALIB.IMG'

roda, file, dx=-15, dy=22, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/N1542049376_1_CALIB.IMG'

roda, file, dx=0, dy=3, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/N1555775942_1_CALIB.IMG'

roda, file, dx=5, dy=-2, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/N1583788619_1_CALIB.IMG'

roda, file, dx=-2, dy=0, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/N1831871631_1_CALIB.IMG'

roda, file, dx=14, dy=0, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/N1546720830_1_CALIB.IMG'

roda, file, dx=-22, dy=0, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/W1467350987_2_CALIB.IMG'

roda, file, dx=10, dy=0, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE

file = 'imagens/W1477738447_1_CALIB.IMG'

roda, file, dx=0, dy=0, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IENCKE, /OENCKE

file = 'prometheus_pandora/N1849967837_1_CALIB.IMG'

roda, file, dx=0, dy=-50, rangex=30, rangey=30, incrx=0.25, incry=0.25, output=matrix, /ARING, /IKEELER, /OKEELER, /IENCKE, /OENCKE, high=0.993, low=0.98

