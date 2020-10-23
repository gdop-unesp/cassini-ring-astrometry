# coding: utf-8
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import gridspec
import numpy as np
import matplotlib.image as mpimg
from mpl_toolkits.mplot3d import Axes3D
from astropy.io import fits
import astropy.units as u
from astropy.visualization import ZScaleInterval, LinearStretch, MinMaxInterval, SqrtStretch, ImageNormalize
import variaveis as var

plt.rcParams['text.latex.preamble']=[r"\usepackage{txfonts}"]

params = {'text.usetex': True,
          'font.size':22,
          'font.family':'txfonts',
          'text.latex.unicode':True,}
sizel=12

fiarin = var.fiarin.T
fiein = var.fiein.T
fieout = var.fieout.T
#fikin = var.fikin.T
#fikout = var.fikout.T
fisats = var.fisats.T
inarin = var.inarin.T
inein = var.inein.T
ineout = var.ineout.T
#inkin = var.inkin.T
#inkout = var.inkout.T
insats = var.insats.T

sats = ['Pan', 'Daphnis', 'Atlas', 'Prometheus', 'Pandora', 'Epimetheus', 'Janus', 'Aegaeon', 'Mimas', 'Methone', 'Anthe', 'Pallene', 'Enceladus', 'Tethys', 'Telesto', 'Calypso', 'Polydeuces', 'Dione', 'Helene', 'Rhea', 'Titan', 'Hyperion', 'Iapetus', 'Phoebe']

an = np.linspace(0,2*np.pi, 100)
r = 3

def color_map(rangex, stepx, dx, rangey, stepy, dy, ticks=5):
    img = fits.getdata('matrix.fits')
    x,y = np.indices(img.shape)

    fig, ax = plt.subplots(figsize=((17.6*u.cm).to(u.imperial.inch).value, (17.6*u.cm).to(u.imperial.inch).value*(9.0/8.0)))
    cax = ax.imshow(img/img.max(), cmap='YlOrRd', origin='lower')
    xt = np.linspace(0,rangex/stepx,ticks)
    xl = xt*stepx - rangex/2.0 + dx
    yt = np.linspace(0,rangey/stepy,ticks)
    yl = yt*stepy - rangey/2.0 + dy
    ax.set_xticks(xt)
    ax.set_xticklabels(xl)
    ax.set_yticks(yt)
    ax.set_yticklabels(yl)
    ax.tick_params(labelsize=sizel)
    ax.set_xlabel('Displacement in X (pixels)', fontsize=sizel)
    ax.set_ylabel('Displacement in Y (pixels)', fontsize=sizel)
    cbar = fig.colorbar(cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=sizel)
    plt.savefig('W1477738447_colorbar.png', dpi=300, bbox_inches='tight')
    plt.close()
    
def img_reg(x1, tx, y1):
    x2 = x1 + tx
    y2 = int(tx*(2/3) + y1)
    print(x1,x2,y1,y2)
    img = fits.getdata('imagem.fits')
    x,y = np.indices(img.shape)
    fig = plt.figure(figsize=((17.6*u.cm).to(u.imperial.inch).value, (17.6*u.cm).to(u.imperial.inch).value*((y2-y1)/(x2-x1))))
    norm = ImageNormalize(img[y1:y2,x1:x2], interval=ZScaleInterval(), stretch=SqrtStretch())
    axf = fig.add_axes([0.0,0.0,1.0,1.0])
    axf.imshow(img[y1:y2,x1:x2], norm=norm, cmap='gray', origin='lower')
    for v in [inarin,inein,ineout]:
        if len(v) > 0:
            k = np.where((v[0] > x1) & (v[0] < x2-1) & (v[1] > y1) & (v[1] < y2-1))
            axf.plot(v[0][k]-x1, v[1][k]-y1, color='red')
    for n,v in enumerate(insats[0]):
        if insats[0][n] > x1 and insats[0][n] < x2-1 and insats[1][n] > y1 and insats[1][n] < y2-1:
            axf.plot(r*np.cos(an)+insats[0][n]-x1,r*np.sin(an)+insats[1][n]-y1, 'r')
    for v in [fiarin,fiein,fieout]:
        if len(v) > 0:
            k = np.where((v[0] > x1) & (v[0] < x2-1) & (v[1] > y1) & (v[1] < y2-1))
            axf.plot(v[0][k]-x1, v[1][k]-y1, color='yellow')
    for n,v in enumerate(fisats[0]):
        if fisats[0][n] > x1 and fisats[0][n] < x2-1 and fisats[1][n] > y1 and fisats[1][n] < y2-1:
            axf.plot(r*np.cos(an)+fisats[0][n]-x1,r*np.sin(an)+fisats[1][n]-y1, 'yellow')
            axf.text(fisats[0][n]-x1+5,fisats[1][n]-y1+5, sats[n], color='yellow', fontsize=sizel)
    plt.savefig('W1477738447.png', dpi=300)
    plt.close()
    
def both(x1, tx, y1, rangex, stepx, dx, rangey, stepy, dy, ticks=5):
    img = fits.getdata('matrix.fits')
    x,y = np.indices(img.shape)
    x2 = x1 + tx
    y2 = int(tx*(2/3) + y1)

    fig = plt.figure(figsize=((17.6*u.cm).to(u.imperial.inch).value, (17.6*u.cm).to(u.imperial.inch).value*((9.0/8.0) + (y2-y1)/(x2-x1))/2.0))
    gs = gridspec.GridSpec(1,2,width_ratios=[1,2])
    #fig, ax = plt.subplots(figsize=((17.6*u.cm).to(u.imperial.inch).value, (17.6*u.cm).to(u.imperial.inch).value*(9.0/8.0)))
    ax0 = plt.subplot(gs[0])
    cax = ax0.imshow(img/img.max(), cmap='YlOrRd', origin='lower')
    xt = np.linspace(0,rangex/stepx,ticks)
    xl = xt*stepx - rangex/2.0 + dx
    yt = np.linspace(0,rangey/stepy,ticks)
    yl = yt*stepy - rangey/2.0 + dy
    ax0.set_xticks(xt)
    ax0.set_xticklabels(xl)
    ax0.set_yticks(yt)
    ax0.set_yticklabels(yl)
    ax0.tick_params(labelsize=sizel)
    ax0.set_xlabel('Displacement in X (pixels)', fontsize=sizel)
    ax0.set_ylabel('Displacement in Y (pixels)', fontsize=sizel)
    #plt.colorbar(cax=cax, orientation='horizontal')
    cbar = plt.colorbar(cax, orientation='horizontal')
    cbar.ax.tick_params(labelsize=sizel)
    
    img = fits.getdata('imagem.fits')
    x,y = np.indices(img.shape)
    #fig = plt.figure(figsize=((17.6*u.cm).to(u.imperial.inch).value, (17.6*u.cm).to(u.imperial.inch).value*((y2-y1)/(x2-x1))))
    ax1 = plt.subplot(gs[1])
    norm = ImageNormalize(img[y1:y2,x1:x2], interval=ZScaleInterval(), stretch=SqrtStretch())
    #axf = ax1.add_axes([0.0,0.0,1.0,1.0])
    ax1.imshow(img[y1:y2,x1:x2], norm=norm, cmap='gray', origin='lower')
    for v in [inarin,inein,ineout]:
        if len(v) > 0:
            k = np.where((v[0] > x1) & (v[0] < x2-1) & (v[1] > y1) & (v[1] < y2-1))
            ax1.plot(v[0][k]-x1, v[1][k]-y1, color='red')
    for n,v in enumerate(insats[0]):
        if insats[0][n] > x1 and insats[0][n] < x2-1 and insats[1][n] > y1 and insats[1][n] < y2-1:
            ax1.plot(r*np.cos(an)+insats[0][n]-x1,r*np.sin(an)+insats[1][n]-y1, 'r')
    for v in [fiarin,fiein,fieout]:
        if len(v) > 0:
            k = np.where((v[0] > x1) & (v[0] < x2-1) & (v[1] > y1) & (v[1] < y2-1))
            ax1.plot(v[0][k]-x1, v[1][k]-y1, color='yellow')
    for n,v in enumerate(fisats[0]):
        if fisats[0][n] > x1 and fisats[0][n] < x2-1 and fisats[1][n] > y1 and fisats[1][n] < y2-1:
            ax1.plot(r*np.cos(an)+fisats[0][n]-x1,r*np.sin(an)+fisats[1][n]-y1, 'yellow')
            ax1.text(fisats[0][n]-x1+5,fisats[1][n]-y1+5, sats[n], color='yellow', fontsize=sizel)
    plt.savefig('teste.png', dpi=300)
    plt.close()
    

#fig = plt.figure()
#ax = fig.gca(projection='3d')
#ax.plot_surface(x,y,img, cmap='coolwarm')
#plt.show()

color_map(30.0, 0.25, 0.0, 30.0, 0.25, 0.0)
img_reg(315,300,540)
#both(1,999,1,30,0.25,22,30,0.25,0)

