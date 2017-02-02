import matplotlib
matplotlib.use('Agg')
from astLib import astWCS, astCoords
import liteMap as lm
import ConfigParser

config = ConfigParser.SafeConfigParser()
config.read("input/general.ini")

savePath = config.get('sims','root_location') + "output/"

N = 38

for i in range(1,N+1):
    print i

    lmap = lm.liteMapFromFits(savePath+"kappa"+str(i).zfill(3)+".fits")
    if i==1:
        mf = lmap.data*0.



    mf += lmap.data


mf = mf / (N-1)

lmap = lm.liteMapFromFits(savePath+"kappa000.fits")
from orphics.tools.output import Plotter

pl = Plotter()
pl.plot2d(mf)
pl.done("mf.png")


pl = Plotter()
pl.plot2d(lmap.data - mf)
pl.done("mfsub.png")

