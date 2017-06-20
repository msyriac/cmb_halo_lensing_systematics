print "Starting imports..."
import matplotlib
matplotlib.use('Agg')
from orphics.analysis.quadraticEstimator import Estimator
import numpy as np
from astLib import astWCS, astCoords
import liteMap as lm
from orphics.tools.output import Plotter
from orphics.tools.stats import binInAnnuli
import sys
import orphics.analysis.flatMaps as fmaps

import pyfits
import ConfigParser

from orphics.tools.cmb import loadTheorySpectraFromCAMB


config = ConfigParser.SafeConfigParser()
config.read("input/general.ini")

fileRoot = config.get('sims','root_location') + config.get('sims','sim_root')

savePath = config.get('sims','root_location') + "output/"

istart = config.getint('sims','istart')
iend = config.getint('sims','iend')

cambRoot = config.get('recon','cambRoot')
TCMB = config.getfloat('recon','TCMB')
theory = loadTheorySpectraFromCAMB(cambRoot,unlensedEqualsLensed=True,useTotal=False,TCMB = TCMB,lpad=9000)

ellkk = np.arange(2,9000,1)
Clkk = theory.gCl("kk",ellkk)    

beamArcmin = config.getfloat('recon','beamArcmin')
noiseT = config.getfloat('recon','noiseT')
noiseP = config.getfloat('recon','noiseP')
cmbellmax = config.getint('recon','cmbellmax')
kellmax = config.getint('recon','kellmax')
tapWid = config.getint('recon','tapWidth')
tapPad = config.getint('recon','tapPad')

avg = 0.
avg2 = 0.

k = 0
for i in range(istart,iend+1):

    froot = fileRoot + str(i).zfill(3)
    imgFile = froot + "_tSZ_MOD.fits"
    hd = pyfits.open(imgFile)
    bigMap = hd[0].data-hd[6].data-hd[7].data
    print bigMap.shape


    
    if k==0:
        wcs = astWCS.WCS(imgFile,extensionName = 0)
        
        lmap = lm.liteMapFromDataAndWCS(bigMap,wcs)
        lmap.info()
        px = np.abs(lmap.x1-lmap.x0)/lmap.Nx*np.pi/180.\
             *np.cos(np.pi/180.*0.5*(lmap.y0+lmap.y1))*180./np.pi*60.
        print px
        lmap.pixScaleX = px/180.*np.pi/60.

        templateMap = lmap.copy()
        lxMap,lyMap,modLMap,thetaMap,lx,ly  = fmaps.getFTAttributesFromLiteMap(templateMap)

        window = fmaps.initializeCosineWindow(templateMap,lenApod=tapWid,pad=tapPad)



        print "Making white noise..."
        nT,nP = fmaps.whiteNoise2D([noiseT,noiseP],beamArcmin,modLMap,TCMB=TCMB)
        fMask = fmaps.fourierMask(lx,ly,modLMap,lmin=2,lmax=cmbellmax)
        fMaskK = fmaps.fourierMask(lx,ly,modLMap,lmin=2,lmax=kellmax)
        qest = Estimator(templateMap,
                         theory,
                         theorySpectraForNorm=None,
                         noiseX2dTEB=[nT,nP,nP],
                         noiseY2dTEB=[nT,nP,nP],
                         fmaskX2dTEB=[fMask]*3,
                         fmaskY2dTEB=[fMask]*3,
                         fmaskKappa=fMaskK,
                         doCurl=False,
                         TOnly=True,
                         halo=True,
                         gradCut=10000,verbose=True)


        # CHECK THAT NORM MATCHES HU/OK
        data2d = qest.AL['TT']
        modLMap = qest.N.modLMap
        bin_edges = np.arange(2,kellmax,10)
        centers, Nlbinned = binInAnnuli(data2d, modLMap, bin_edges)

        huFile = '/astro/u/msyriac/repos/cmb-lensing-projections/data/NoiseCurvesKK/hu_tt.csv'
        huell,hunl = np.loadtxt(huFile,unpack=True,delimiter=',')

        pl = Plotter(scaleY='log',scaleX='log')
        pl.add(ellkk,4.*Clkk/2./np.pi)
        pl.add(centers,4.*Nlbinned/2./np.pi)#,ls="none",marker="o")
        pl.add(huell,hunl,ls='--')#,ls="none",marker="o")
        pl.done("testbin.png")

    passMap = bigMap[:,:]*window[:,:]
    passMap = passMap - passMap.mean()

    if k==0:
        pl = Plotter()
        pl.plot2d(passMap)
        pl.done("bigmap.png")


    print "Reconstructing" , i , " ..."
    qest.updateTEB_X(passMap.astype(float)/2.7255e6)
    qest.updateTEB_Y(passMap.astype(float)/2.7255e6)

    kappa = qest.getKappa('TT')
    kappa = kappa - kappa.mean()

    if k==0:
        pl = Plotter()
        pl.plot2d(kappa.real)
        pl.done("kappa.png")

        #meanfield = kappa.real.copy()*0.


    w2 = np.mean(window**2.)
    templateMap.data = kappa.real/w2
    templateMap.writeFits(savePath+"kappa"+str(i).zfill(3)+".fits",overWrite=True)

    lmap1 = lm.liteMapFromFits(savePath+"inputKappa"+str(i).zfill(3)+".fits")

    if k==0:
        pl = Plotter()
        pl.plot2d(lmap1.data)
        pl.done("inkappa.png")



    px = np.abs(lmap1.x1-lmap1.x0)/lmap1.Nx*np.pi/180.\
         *np.cos(np.pi/180.*0.5*(lmap1.y0+lmap1.y1))*180./np.pi*60.
    print px
    lmap1.pixScaleX = px/180.*np.pi/60.

    lmap1.data = lmap1.data*window
    w4 = np.mean(window**4.)

    print "crossing with input"
    import fftTools as ft
    import orphics.tools.stats as stats
    p2d = ft.powerFromLiteMap(templateMap,lmap1,applySlepianTaper=False)
    bin_edges = np.arange(100,2000,50)
    centers, means = stats.binInAnnuli(p2d.powerMap, p2d.modLMap, bin_edges)

    avg += means



    pl = Plotter()
    pl.add(centers,centers*avg/k/w4)
    pl.add(ellkk,ellkk*Clkk,lw=2)
    pl._ax.set_xlim(100,2000)
    pl.done("crosspower.png")

    p2d = ft.powerFromLiteMap(lmap1,applySlepianTaper=False)
    centers, means = stats.binInAnnuli(p2d.powerMap, p2d.modLMap, bin_edges)

    avg2 += means

    w = 1.#np.sqrt(np.mean(window))

    pl = Plotter()
    pl.add(centers,centers*avg2/k/w4,ls="none",marker="o",markersize=10)
    pl.add(ellkk,ellkk*Clkk,lw=2)
    pl._ax.set_xlim(100,2000)
    pl.done("autopower.png")

    k+=1


