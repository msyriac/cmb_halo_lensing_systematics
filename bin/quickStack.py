import pyfits
import ConfigParser
import numpy as np
import sys


def main(argv):

    config = ConfigParser.SafeConfigParser()
    config.read("input/general.ini")

    fileRoot = config.get('sims','root_location') + config.get('sims','sim_root')

    istart = config.getint('sims','istart')
    iend = config.getint('sims','iend')


    widthArcmin = config.getfloat('cutouts','widthArcmin')
    pixScaleArcmin = config.getfloat('sims','pixScaleArcmin')

    npixx = int(widthArcmin/pixScaleArcmin/2.)
    npixy = int(widthArcmin/pixScaleArcmin/2.)

    #print npixx, npixy

    stack = 0.

    tot = 0
    skip = 0
    for i in range(istart,iend+1):

        froot = fileRoot + str(i).zfill(3)

        catFile =  froot + ".txt"
        imgFile = froot + "_tSZ_MOD.fits"


        xcens, ycens, m200s = np.loadtxt(catFile,unpack=True,usecols=[1,2,6])

        selection = m200s>1.e14

        print "Loading fits file ", i, "with ", len(xcens[selection]), " clusters ..."
        hd = pyfits.open(imgFile)

        bigMap = hd[0].data
        bigY,bigX = bigMap.shape
        print "Map area ", bigX*bigY*pixScaleArcmin*pixScaleArcmin/60./60. , " sq. deg."


        for xcen,ycen,m200 in zip(xcens[selection],ycens[selection],m200s[selection]):
            tot += 1
            xcen = int(xcen)
            ycen = int(ycen)

            if ycen-npixy<0 or xcen-npixx<0 or ycen+npixy>(bigY-1) or xcen+npixx>(bigX-1):
                skip += 1
                continue

            cutOut = bigMap[(ycen-npixy):(ycen+npixy) , (xcen-npixx):(xcen+npixx)]
        

            stack += cutOut


    print "Percentage near bounds = ", skip*100./tot , " %"
    from orphics.tools.output import Plotter
    
    pl = Plotter()
    pl.plot2d(stack)
    pl.done("plots/stack.png")




if (__name__ == "__main__"):
    main(sys.argv[1:])

