from astropy.io import fits
import numpy as np
import glob


folderSub = '/home/oan1803/Documents/MOAT/Mancha1/'

files = sorted(glob.glob(folderSub+'*.fits'))

imagen0 = fits.getdata(files[312])

media0 = np.mean(imagen0)
desvi0 = np.std(imagen0)
diferencia = media0 - 3*desvi0


for i in range(len(files)):
	
	imageni = fits.getdata(files[i])
	media = np.mean(imageni)
	# ~ desvi = np.std(imageni)
	
	
	if media < diferencia: 
		print(i-1)
		break
