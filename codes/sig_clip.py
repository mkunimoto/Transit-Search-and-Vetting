# Sigma-clip a light curve in the positive flux direction

import sys
import numpy as np

data = np.loadtxt(sys.argv[1])
sig = float(sys.argv[2])

time = data[:,0]
flux = data[:,1]
err = data[:,2]

iter = 5

i = 0
while i < iter:
	i += 1 
	MAD = np.median(abs(flux - np.median(flux)))
	std = 1.48*MAD
	time = time[(flux < sig*std)]
	err = err[(flux < sig*std)]
	flux = flux[(flux < sig*std)]

new_data = np.column_stack((time,flux,err))

np.savetxt("sig_clip.dat",new_data)
