# -*- coding: utf-8 -*-
"""
This example demonstrates the use of ImageView, which is a high-level widget for 
displaying and analyzing 2D and 3D data. ImageView provides:

  1. A zoomable region (ViewBox) for displaying the image
  2. A combination histogram and gradient editor (HistogramLUTItem) for
     controlling the visual appearance of the image
  3. A timeline for selecting the currently displayed frame (for 3D data only).
  4. Tools for very basic analysis of image data (see ROI and Norm buttons)

"""
## Add path to library (just for examples; you do not need this)
#~ import initExample
import matplotlib
matplotlib.use('Qt5Agg')
import numpy as np
from pyqtgraph.Qt import QtCore, QtGui
import pyqtgraph as pg
from astropy.io import fits
import matplotlib.pyplot as plt
from cmapToColormap import *
from scipy.io import readsav

# Interpret image data as row-major instead of col-major
pg.setConfigOptions(imageAxisOrder='row-major')

app = QtGui.QApplication([])

## Create window with ImageView widget
win = QtGui.QMainWindow()
win.resize(800,800)
imv = pg.ImageView(view=pg.PlotItem())
win.setCentralWidget(imv)
win.show()
win.setWindowTitle('pyqtgraph example: ImageView')

## Create random 3D data set with noisy signals
#~ img = pg.gaussianFilter(np.random.normal(size=(200, 200)), (5, 5)) * 20 + 100
#~ img = img[np.newaxis,:,:]
#~ decay = np.exp(-np.linspace(0,0.3,100))[:,np.newaxis,np.newaxis]
#~ data = np.random.normal(size=(100, 200, 200))
#~ data += img * decay
#~ data += 2

#~ ## Add time-varying signal
#~ sig = np.zeros(data.shape[0])
#~ sig[30:] += np.exp(-np.linspace(1,10, 70))
#~ sig[40:] += np.exp(-np.linspace(1,10, 60))
#~ sig[70:] += np.exp(-np.linspace(1,10, 30))

#~ sig = sig[:,np.newaxis,np.newaxis] * 3
#~ data[:,50:60,30:40] += sig

# ~ data2 = fits.getdata("/Users/joseivan/Desktop/KSO/cube_algine.fits")
# ~ data2 = np.swapaxes(np.swapaxes(data2, 1, 2), 0, 1)
print('Loading data')
# ~ data2 = fits.getdata("/Volumes/DI/HIFI_DATA/Cube/cube_hifi_20170928_alinea_2_filtered_speckle.fts")
data2 = fits.getdata('/Volumes/VALERIA/AUSUS/Documents/MOAT/Cubos/cubealign3hours.fits')
# ~ dict_V = readsav('/Users/joseivan/Documents/DataGREGOR/28sep17_V.sav')
# ~ vv = dict_V['vv']
# ~ del(dict_V)
# ~ data2 = np.swapaxes(vv,1,0)

# ~ data2 = np.swapaxes(data2,0,3)
# ~ data2 = (data2[0:365,2,:,:]).astype(np.float16)
print(data2.dtype,data2.flags.owndata,data2.shape)
#data2 = np.lib.format.open_memmap('/Users/joseivan/Desktop/KSO/cube_aligne.npy', dtype=np.int16, mode='r')
print('Normalizing data')
#data2 = data2[:,::-1,:]/data2.max()
## Display the data and assign each frame a time value from 1.0 to 3.0
print('creating plot...')
# ~ imv.setImage(data2[:,::-1,:], xvals=np.linspace(0, 360,data2.shape[0]))
imv.setImage(data2[:,:,:], xvals=np.linspace(0, data2.shape[0],data2.shape[0]))

## Set a custom color map
#~ colors = [
    #~ (0, 0, 0),
    #~ (45, 5, 61),
    #~ (84, 42, 55),
    #~ (150, 87, 60),
    #~ (208, 171, 141),
    #~ (255, 255, 255)
#~ ]
#~ cmap = pg.ColorMap(pos=np.linspace(0.0, 1.0, 6), color=colors)
pos, rgba_colors = zip(*cmapToColormap(plt.get_cmap('gray')))
        # Set the colormap
cmap =  pyqtgraph.ColorMap(pos, rgba_colors)
imv.setColorMap(cmap)

## Start Qt event loop unless running in interactive mode.
if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
