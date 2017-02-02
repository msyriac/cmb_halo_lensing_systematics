import matplotlib
matplotlib.use('Agg')
import ConfigParser
import numpy as np
import matplotlib.pyplot as plt
import glob
class BattagliaSims(object):

    def __init__(self,rootPath="/astro/astronfs01/workarea/msyriac/clusterSims/Battaglia/"):

        self.root = rootPath

        

    def snapToZ(self,snapNum):

    def getTSZ(self):
        PIX = 2048 # all files are 2048 x 2048
        filelist = glob.glob(self.root+"")
        i = 0
        for filen in filelist:
            with open(filen, 'rb') as fd:
                temp = np.fromfile(file=fd, dtype=np.float32)

            #reshape array into 2D array
            map = np.reshape(temp,(PIX,PIX))
            temp = 0

            #Check map values to make sure they aren't crazy
            print np.max(map), np.min(map)
            #Check the image
            zoom = PIX/2 -128
            plt.imshow((map[zoom:-zoom,zoom:-zoom]),interpolation='nearest')
            plt.savefig("plots/plot"+str(i)+".png")
            i+=1



b  =  BattagliaSims()
b.getTSZ()
