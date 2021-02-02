# sdo_crop
import time
import numpy as np
import glob
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.coordinates import Helioprojective
from sunpy.physics.differential_rotation import differential_rotate
from sunpy.time import TimeRange
import sunpy.map as smap



print('Starting the crop code')
t0 = time.time()
folder = '/media/hypnus/VALERIA/DAVELCT/CompleteField/'
files = sorted(glob.glob(folder+'*.fits'))

indx_0 = 0
indx_f = indx_0+len(files)

indx_c = 8 #Se escoge el dato de la mitad para alinearlo con el 


xc = -113*u.arcsec
yc = -235.1*u.arcsec
length = 90*u.arcsec

xmap_c = smap.Map(files[indx_c])
time_c = xmap_c.meta['date-obs']

folder_test = '/media/hypnus/VALERIA/DAVELCT/Submapfield/'

for i in range(indx_0,indx_f):
	print('Reducing frame --> ',i)
	
	xmap = smap.Map(files[i])
	time_ = xmap.meta['date-obs']
	duration = TimeRange(time_,time_c).seconds
		
	
	rmap = xmap.rotate()
	
	if (indx_c - i) >= 0:
		dmap = differential_rotate(rmap, time=duration)
		
	elif (indx_c - i) < 0:
		dmap = differential_rotate(rmap, time=-duration)
	

	bottom_left = SkyCoord(xc - length, yc - length,
							   frame=dmap.coordinate_frame)
	top_right = SkyCoord(xc + length, yc + length,
							 frame=dmap.coordinate_frame)
	
	submap = dmap.submap(bottom_left, top_right)
	
	submap.save(folder_test+'submap_'+str(i).zfill(4)+'.fits',overwrite=True)
print('Time elapsed creating full submaps --> ',time.time()-t0)


