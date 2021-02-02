#!/usr/bin/env python
# -*- coding: utf8 -*-
#EL codigo hace un filtro pasa bajos, quita todo lo que esta por debajo de cierta velocidad.
# ~ import pyfits as fits
from astropy.io import fits
import numpy as np
import os
import time
from bignfft import *
import shutil
import glob
import sys
import pickle


print("Iniciando codigo")

def bigsonic(cube,first,last,bxdim,bydim):
	if not os.path.exists("/home/oan1803/Documents/MOAT/work"): #Revisa si existe el folder y si no, lo crea (CAMBIA)
		os.makedirs("/home/oan1803/Documents/MOAT/work") #Crea los folders (CAMBIA)
	if not os.path.exists("/home/oan1803/Documents/MOAT/filter"): #(CAMBIA)
		os.makedirs("/home/oan1803/Documents/MOAT/filter") #Crea los folders (CAMBIA)
		
	str0 = "/home/oan1803/Documents/MOAT/work/"     #(CAMBIA)
	strf = "/home/oan1803/Documents/MOAT/filter/"   #(CAMBIA)
	
	
	
	
	
	# lineas editadas para que funcione con el cubo generado por Sunpy
	
	print("Cargando Cubo")
	
	# ====================================================
	# Working with several images created by SunPy
	
	#~ folderOut = "/Volumes/DataIvanOAN/StartingAGAIN/DataPaper/DataSubFields/SubImaIc_cor/"
	#~ subs = sorted(glob.glob(folderOut+"*.fits"))
	
	#~ cube = smap.Map(subs,cube=True)
	
	#~ maps = cube.maps
	
	#~ zdim = len(maps)
	
	#~ dimx = maps[0].data.shape[1]
	#~ dimy = maps[0].data.shape[0]
	# ===================================================
	
	
	# ====================================================
	# Working with cube data made by pyfits or astropy
	
	# ~ cubename = input("Please, put the full path of the cube:")
	# ~ cube =fits.getdata(cubename)
	
	#cube= "/home/oan1803/Documents/MOAT/Cubos/cubito.fits"
	dimx = cube.shape[2]
	dimy = cube.shape[1]
	# ===================================================
	
	
	
	xdim=dimx
	ydim=dimy
	x_anf=0
	y_anf=0
	scale = 0.504 #arcsec/pixel ; ** change manually **
	
#-----------------------------------------------------------------

	if xdim%2 != 0: xdim = xdim - 1	
	if ydim%2 != 0: ydim = ydim - 1	

	if (last - first + 1)%2 != 0: last = last - 1	

	tdim = last - first + 1

	print("-------")
	#Mean time separation between images [s] ; ** change manually **
	t_step = 720
	
	# Maximum phase velocity [km/s]
	v_ph = 4.
	
	ap = 0
	
	cut = 0
	
	if (cut > 2) or (cut < 0):
		raise ValueError('The cut values alowed are 0, 1, and 2')
	if (ap != 0) or (cut != 0):
		perct = 10
		smooth_t = int(tdim*perct/100) # width of the edge in t-dimension

	t1 = time.time()
	
	#-----------------------------------------
	
	print(tdim, " images to be filtered")
	
	# subsonic construction
	
	#Definin' unit in the Fourier domain
	
	print("Spatial resolution ->", scale, " arcsec/pixel")
	
	print("---")
	
	kx_step = 1./(scale*725.*xdim)
	ky_step = 1./(scale*725.*ydim)
	w_step = 1./(t_step*tdim)
	
	
	# Prepare de filter
	
	nx = int(xdim/2) + 1
	ny = int(ydim/2) + 1
	nt = int(tdim/2) + 1
	filter_mask = np.zeros([ny,nx,nt],dtype=float)
	
	
	if (cut == 0) or (cut == 1):
		print("Now calculatin' filter...")
		for j in range(ny):
			for i in range(nx):
				k_by_v = np.sqrt((i*kx_step)**2+(j*ky_step)**2)*v_ph
				for k in range(1,nt):
					if k*w_step <= k_by_v:
						filter_mask[j,i,k]=1.
	
	if cut == 1:
		for i in range(nx):
			for j in range(ny):
				for k in range(1,nt):
					if filter_mask[j,i,k] != 1:
						trans = perct*k/100
						trans2 = int(trans/2)
						trans = 2*trans
						if trans2 >= 1:
							n=0
							for kk in range(int(k-trans2),int(k+trans2)):
								if kk <= nt-1:
									filter_mask[j,i,kk] = 0.5+0.5*np.cos(np.pi*n/trans)
									n = n+1
								else:
									continue
	

	
	for i in range(nt):
		filter_slice = np.zeros([ydim,xdim])
		if i == 0:
			filter_slice = filter_slice + 1.
		else:
			filter_slice[0:ny,0:nx]=filter_mask[:,:,i]
			filter_slice[ny:ydim,:]=filter_slice[1:ny-1,:][::-1,:]
			filter_slice[:,nx:xdim]=filter_slice[:,1:nx-1][:,::-1]
		dcn = str(i).zfill(4)
		pickle.dump( filter_slice, open( strf+dcn+".p", "wb" ) )
		
		if i != 0:
			dcn = str(tdim-i+first).zfill(4)
			pickle.dump( filter_slice, open( strf+dcn+".p", "wb" ) )

	print("The filter was written to "+strf)
	
	del(filter_mask)
	del(filter_slice)
	
	
	if ap != 0:
		tmask = np.ones(tdim)
		for i in range(smooth_t):
			tmask[i] = (1-np.cos(np.pi*i/smooth_t))/2
		tmask[tdim-smooth_t:tdim] = (tmask[1:smooth_t+1])[::-1]
		print('Computing the mean for the cube (it could take several minutes):')
		av = 0.
		for n in range(first,last+1):
			ima = cube[n,:,:]
			ima = ima[y_anf:y_anf+ydim,x_anf:x_anf+xdim]
			av = av+np.mean(ima)/tdim
	
	
	# Loop of reading, optional apodization and writing images
	for n in range(first,last+1):
		dcn = str(n).zfill(4)
		print("Reding images -->", "cube["+str(n)+"]")
		
		# Editado para trabajar con sunpy
		
		#~ ima = maps[n].data
		ima = cube[n,:,:]
		ima = ima.astype('float64')
		
		# apodization
		if ap != 0:
			ima = ima-av
			ima = ima*tmask[n-first]
			ima = ima + av
			del(tmask)
		

		pickle.dump( ima, open( str0+"apo"+dcn+".p", "wb" ) )
		
		print("apodized images were written to "+str0)
	
	del(ima)
	
	
	# Direct FFT
	
	print(50*"=")
	print("Calling bignfft")
	print(50*"-")
	
	
	bignfft(-1, first,last,dimx,dimy,bxdim,bydim)
	
	
	
	# Multiplying transformed images by filter images
	
	
	for n in range(first,last+1):
		dcn = str(n).zfill(4)
		ima = pickle.load( open( str0+"fft"+dcn+".p", "rb" ) )
		filter = pickle.load( open( strf+dcn+".p", "rb" ) )
		
		print("multiplying by the filter ",str0+"fft"+dcn)
		
		ima = ima*filter
		pickle.dump( ima, open( str0+"fft"+dcn+".p", "wb" ) )

	del(ima)
	del(filter)
	
	#Inverse FFT
	print(50*"=")
	print("Calling bignfft")
	print(50*"-")
	
	bignfft(1, first,last,dimx,dimy,bxdim,bydim)

	
	# Saving results
	
	ima = np.zeros([ydim,xdim],dtype="complex64")
	
	cube_new = np.zeros([tdim,ydim,xdim],dtype='float64')
	
	print((ima.real).shape)
	print(cube_new.shape)
	
	
	for n in range(first,last+1):
		dcn = str(n).zfill(4)
		ima = pickle.load( open( str0+"F"+dcn+".p", "rb" ) )
		im = (ima.real).astype('float64')
		cube_new[n,:,:] = im
	
	# ~ hdu = fits.PrimaryHDU(cube_new)
	# ~ hdulist = fits.HDUList([hdu])
	# ~ hdulist.writeto(str2+'cubeQuietSun_Mag_filtered_20110411.fits',overwrite=True)
	
	print("---")
	print("Total elapsed time from begining = ", np.round(time.time()-t1,2))
	print(" ")
	print("Erasing directories work and filter")
	
	shutil.rmtree(str0)
	shutil.rmtree(strf)
	return cube_new

#bigsonic(0,2080,int(296/15),int(296/15))
		

	
		
	
		
	
	
	
	
	
	
