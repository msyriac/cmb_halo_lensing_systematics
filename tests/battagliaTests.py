import matplotlib
matplotlib.use('Agg')
from ConfigParser import SafeConfigParser 
import numpy as np
from szlib.sims import BattagliaSims
from orphics.tools.io import Plotter,dictFromSection,listFromConfig


# === COSMOLOGY ===
iniFile = "../SZ_filter/input/params.ini"
Config = SafeConfigParser()
Config.optionxform=str
Config.read(iniFile)
lmax = 5000
constDict = dictFromSection(Config,'constants')


b  =  BattagliaSims(constDict)
bsims = b.mapReader()

for ret in bsims:
    z,mapList = ret
    dm,star,gas,sz = mapList
    # print dm.mean()
    # print star.mean()
    # print gas.mean()
    # print sz.mean()
