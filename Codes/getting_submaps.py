import glob
import numpy as np
from mpi4py import MPI
import sunpy.map as smap
import astropy.units as u
from astropy.coordinates import SkyCoord
from sunpy.physics.differential_rotation import solar_rotate_coordinate
from sunpy.time import parse_time
import sys


print('Starting the reduction..')

folder_in = '/home/oan1803/Documents/MOAT/Mancha_Raw/Mancha_estrella/'
folder_out = '/home/oan1803/Documents/MOAT/Mancha1/'

files = sorted(glob.glob(folder_in+'*.fits'))
ind_c = int(np.median(np.arange(len(files))))
length = 140*u.arcsec

# ~ print('submap'+files[ind_c][60:])
# ~ sys.exit()

print('Creating central map')
	
map_c = smap.Map(files[ind_c])
rot_c = map_c.rotate()

xc , yc = 216*u.arcsec, 95*u.arcsec

time_c = rot_c.meta['date-obs']
start_coord = SkyCoord(xc, yc, frame=rot_c.coordinate_frame)

bottom_left_c = SkyCoord(start_coord.Tx-length/2,start_coord.Ty-length/2,frame=rot_c.coordinate_frame)
top_right_c = SkyCoord(start_coord.Tx+length/2,start_coord.Ty+length/2,frame=rot_c.coordinate_frame)
submap_c = rot_c.submap(bottom_left_c, top_right_c)
submap_c.save(folder_out+'submap'+files[ind_c][59:])
print('Central submap created...')
	

for i in range(len(files)):
	print('Reducing file {0}'.format(i))
	
	if i != ind_c:
			
		map_i = smap.Map(files[i])
		rot_i = map_i.rotate()
		time_i = rot_i.meta['date-obs']
		diff_coord = solar_rotate_coordinate(start_coord, time=parse_time(time_i))
		
		bottom_left = SkyCoord(diff_coord.Tx-length/2,diff_coord.Ty-length/2,frame=rot_i.coordinate_frame)
		top_right = SkyCoord(diff_coord.Tx+length/2,diff_coord.Ty+length/2,frame=rot_i.coordinate_frame)
		submap_i = rot_i.submap(bottom_left, top_right)
		
		submap_i.save(folder_out+'submap'+files[i][59:])


























