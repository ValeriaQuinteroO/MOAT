 #!/usr/bin/env python
# -*- coding: utf8 -*-


from astropy.io import fits
import time
import glob
import numpy as np


# Making cubes...

print("Creating a cube data with correct dimensions and a pseudo tracking")


folderSub = "/home/oan1803/Documents/MOAT/Mancha1/"
folderCubes = "/home/oan1803/Documents/MOAT/Cubos/"

#Reading submaps
Submaps = sorted(glob.glob(folderSub+"*.fits"))

t3 = time.time()

a = 0
b = 838



cube = np.zeros([b-a, 277, 277],dtype=float)



for k in range(b-a):
	smap = fits.getdata(Submaps[k+a])
	smap=smap[0:277,0:277]
	# ~ f=smap.shape
	# ~ if f[1]>278:
		# ~ print(k,smap.shape)
	print("Adding submap --> ",k+1)
	cube[k,:,:] = smap

# ~ smap = fits.getdata(Submaps[4])	
# ~ print(smap.shape)

# ~ fr=smap[0]

# ~ print(fr[0:2])

# ~ print(smap[0:2,0:2])

hdu = fits.PrimaryHDU(cube)
hdulist = fits.HDUList([hdu])
hdu.writeto(folderCubes+'cubito.fits')

print("Time writing cube -->", time.time()-t3)
