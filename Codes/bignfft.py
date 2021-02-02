#!/usr/bin/env python
# -*- coding: utf8 -*-

#~ import pyfits as fits
#~ from astropy.io import fits
import numpy as np
import os
import time
import sys
from scipy.fftpack import fft,ifft,fft2,ifft2
import pickle


def bignfft(pm, first,last,dimx,dimy,bxdim,bydim):
	str0 = "/home/oan1803/Documents/MOAT/work/"   #(CAMBIA)
	xdim = dimx
	ydim = dimy
	x_inf=0
	y_inf = 0
	if xdim%2 != 0: xdim = xdim - 1	
	if ydim%2 != 0: ydim = ydim - 1

	# --------------------------------------------------------
	
	# Reads separate images and computes their  2d-Fourier transform
	for i in range(first,last+1):
		dcn = str(i).zfill(4)
		
		if pm == -1:
			im = pickle.load( open( str0+"apo"+dcn+".p", "rb" ) )
			fim = fft2(im[y_inf:y_inf+ydim,x_inf:x_inf+xdim])/(dimx*dimy)
			# ~ fim = fft2(im[y_inf:y_inf+ydim,x_inf:x_inf+xdim])/(dimx*dimy)
			# ~ fim = fim.astype(np.complex64)
			print("... now writing its 2d-Fourier transform in "+ str0)
			pickle.dump( fim, open( str0+"fft"+dcn+".p", "wb" ) )
			
		elif pm == 1:
			im = pickle.load( open( str0+"fft"+dcn+".p", "rb" ) )
			fim = ifft2(im[y_inf:y_inf+ydim,x_inf:x_inf+xdim])*(dimx*dimy)
			# ~ fim = ifft2(im[y_inf:y_inf+ydim,x_inf:x_inf+xdim])*(dimx*dimy)
			print("... now writing its 2d-Fourier transform in "+ str0)
			pickle.dump( fim, open( str0+"F"+dcn+".p", "wb" ) )

	del(fim)
	del(im)
	
	
	#Creates 3d subarrays (nu_x,nu_y,t) with dimensions in spatial freq. small
	#enough to enter the 3d subarray in the cpu, and computes the fft of
	#1d temporal arrays corresponding to separated x, y pixels
	
	a = xdim % bxdim
	b = ydim % bydim
	n = int(xdim / bxdim)
	m = int(ydim / bydim)
	
	
	if a != 0:
		dix = np.zeros([n+1],dtype=int)+bxdim
		dix[n]=a
	else:
		dix = np.zeros([n],dtype=int)+bxdim
	
	if b != 0:
		diy = np.zeros([m+1],dtype=int)+bydim
		diy[m] = b
	else:
		diy = np.zeros([m],dtype=int)+bydim
	
	nelx = len(dix)
	nely = len(diy)
	print("Number of subarrays-->", str(nelx*nely))
	
	
	# loop over each sub-box
	
	fim = np.zeros([ydim,xdim],dtype="complex64")
	for jbox in range(nely):
		for ibox in range(nelx):
			box3d = np.zeros([last-first+1,diy[jbox],dix[ibox]],dtype=np.complex64)
			print("Now reading in the Fourier domain the sub-array-->",ibox,jbox)
			for i in range(first,last+1):
				dcn = str(i).zfill(4)
				if pm == -1:
					fim = pickle.load( open( str0+"fft"+dcn+".p", "rb" ) )
				if pm == 1:
					fim = pickle.load( open( str0+"F"+dcn+".p", "rb" ) )
				box3d[i-first,:,:]=fim[jbox*bydim:jbox*bydim+diy[jbox],
				                       ibox*bxdim:ibox*bxdim+dix[ibox]]
			for y in range(diy[jbox]):
				for x in range(dix[ibox]):
					tt = box3d[:,y,x]
					if pm == -1:
						box3d[:,y,x]=fft(tt)/(last-first+1)
						# ~ box3d[:,y,x]=fft(tt)/(last-first+1)
					if pm == 1:
						box3d[:,y,x]=ifft(tt)*(last-first+1)
						# ~ box3d[:,y,x]=ifft(tt)*(last-first+1)
			
			print("... writin' its t-Fourier transform "+str0)
			
			for i in range(first,last+1):
				dcn = str(i).zfill(4)
				if pm == -1:
					fim = pickle.load( open( str0+"fft"+dcn+".p", "rb" ) )
				if pm == 1:
					fim = pickle.load( open( str0+"F"+dcn+".p", "rb" ) )
				
				fim[jbox*bydim:jbox*bydim+diy[jbox],
				                       ibox*bxdim:ibox*bxdim+dix[ibox]] = box3d[i-first,:,:]
				

				
				if pm == -1:
					pickle.dump( fim, open( str0+"fft"+dcn+".p", "wb" ) )

				if pm == 1:
					pickle.dump( fim, open( str0+"F"+dcn+".p", "wb" ) )
				
	print(" ")
	print("=== END of FFT ===")
		
				
		
			
		
		
	
	
	
