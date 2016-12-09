import matplotlib
matplotlib.use('Agg')
import pyfits
import ConfigParser
import numpy as np
import sys

from astLib import astWCS, astCoords
import liteMap as lm

from scipy.fftpack import fft2,ifft2
def main(argv):


    # Read some parameters from the ini file
    config = ConfigParser.SafeConfigParser()
    config.read("input/general.ini")
    fileRoot = config.get('sims','root_location') + config.get('sims','sim_root')
    savePath = config.get('sims','root_location') + "output/"
    istart = config.getint('sims','istart')
    iend = config.getint('sims','iend')
    widthArcmin = config.getfloat('cutouts','widthArcmin')
    pixScaleArcmin = config.getfloat('sims','pixScaleArcmin')

    
    npixx = int(widthArcmin/pixScaleArcmin/2.)
    npixy = int(widthArcmin/pixScaleArcmin/2.)


    stackTot = 0.  # stack on lensed cmb + sz
    stackSZ = 0.   # stack on sz
    stackDx = 0.   # stack on deflection_x
    stackDy = 0.   # stack on deflection_y
    stackKappa = 0.# stack on kappa = div.Deflection / 2
    stackAlt = 0.  # stack on kappa = div.Deflection / 2 calculated with FFTs
    stackEst = 0.  # stack on kappa reconstruction

    totnum = 0  # total number of stamps read
    skip = 0    # number skipped because stamp falls outside edges
    for i in range(istart,iend+1):

        froot = fileRoot + str(i).zfill(3)

        # read catalog of clusters and make a selection on mass
        catFile =  froot + ".txt"
        xcens, ycens, m200s = np.loadtxt(catFile,unpack=True,usecols=[1,2,6])
        selection = m200s>1.e14


        # open cmb map
        imgFile = froot + "_tSZ_MOD.fits"
        print "Loading fits file ", i, "with ", len(xcens[selection]), " clusters ..."
        hd = pyfits.open(imgFile)


        # open kappa reconstruction and apply a dumb correction to the pixel scale
        lmap = lm.liteMapFromFits(savePath+"kappa"+str(i).zfill(3)+".fits")
        px = np.abs(lmap.x1-lmap.x0)/lmap.Nx*np.pi/180.\
             *np.cos(np.pi/180.*0.5*(lmap.y0+lmap.y1))*180./np.pi*60.
        lmap.pixScaleX = px/180.*np.pi/60.
        print "pixel scale arcminutes X , " , lmap.pixScaleX * 180.*60./np.pi
        print "pixel scale arcminutes X , " , lmap.pixScaleY * 180.*60./np.pi



        tot = hd[0].data
        sz = hd[6].data + hd[7].data
        unl = hd[3].data

        dx = hd[8].data*np.pi/180./60.
        dy = hd[9].data*np.pi/180./60.
        kappa = 0.5*(np.gradient(dx,axis=0)+np.gradient(dy,axis=1))
        diffmap = tot - sz - unl
        
        #alt kappa
        import orphics.analysis.flatMaps as fmaps
        lxMap,lyMap,modLMap,thetaMap,lx,ly = fmaps.getFTAttributesFromLiteMap(lmap)
        win = 1.#fmaps.initializeCosineWindow(lmap,lenApod=100,pad=10)
        Ny = lmap.Ny
        kgradx = lx*fft2(dy*win)*1j
        kgrady = ly.reshape((Ny,1))*fft2(dx*win)*1j
        altkappa = 0.5*ifft2(kgradx+kgrady).real
        saveMap = lmap.copy()
        saveMap.data = altkappa
        saveMap.writeFits(savePath+"inputKappa"+str(i).zfill(3)+".fits",overWrite=True)

        if i==0:
            from orphics.tools.output import Plotter
            pl = Plotter()
            pl.plot2d(diffmap*win)
            pl.done("diff.png")
            #sys.exit()

        bigY,bigX = tot.shape
        print "Map area ", bigX*bigY*pixScaleArcmin*pixScaleArcmin/60./60. , " sq. deg."

        doRandom = False
        
        if doRandom:
            xlbound = xcens[selection].min()
            xrbound = xcens[selection].max()
            ylbound = ycens[selection].min()
            yrbound = ycens[selection].max()



        for xcen,ycen,m200 in zip(xcens[selection],ycens[selection],m200s[selection]):
            totnum += 1
            if doRandom:
                xcen = np.random.randint(xlbound,xrbound)
                ycen = np.random.randint(ylbound,yrbound)
            else:
                xcen = int(xcen)
                ycen = int(ycen)

            if ycen-npixy<0 or xcen-npixx<0 or ycen+npixy>(bigY-1) or xcen+npixx>(bigX-1):
                skip += 1
                continue

            cutOutTot = tot[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
            cutOutSZ = sz[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
            cutOutDx = dx[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
            cutOutDy = dy[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
            cutOutKappa = kappa[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
            cutOutAlt = altkappa[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
            cutOutEst = lmap.data[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]

            stackTot += cutOutTot
            stackSZ += cutOutSZ
            stackDx += cutOutDx
            stackDy += cutOutDy
            stackKappa += cutOutKappa
            stackAlt += cutOutAlt
            stackEst += cutOutEst


    print "Percentage near bounds = ", skip*100./totnum , " %"
    from orphics.tools.output import Plotter

    if doRandom:
        rs = "R"
    else:
        rs = ""

    N = totnum - skip

    stackD = np.sqrt(stackDx**2.+stackDy**2.)
    pl = Plotter()
    pl.plot2d(stackD/N)
    pl.done("plots/stackD"+rs+".png")



    pl = Plotter()
    pl.plot2d(stackKappa/N)
    pl.done("plots/stackK"+rs+".png")

    pl = Plotter()
    pl.plot2d(stackAlt/N)
    pl.done("plots/stackA"+rs+".png")

    pl = Plotter()
    pl.plot2d(stackTot/N)
    pl.done("plots/stackTot"+rs+".png")


    pl = Plotter()
    pl.plot2d(stackSZ/N)
    pl.done("plots/stackSZ"+rs+".png")



    pl = Plotter()
    pl.plot2d(stackEst/N)
    pl.done("plots/stackE"+rs+".png")



if (__name__ == "__main__"):
    main(sys.argv[1:])

