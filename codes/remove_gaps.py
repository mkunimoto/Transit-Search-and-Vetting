# Remove data near gaps in a light curve

import sys
import numpy as np

data = np.loadtxt(sys.argv[1])
time = data[:,0]
flux = data[:,1]
err = data[:,2]

gap_size = float(sys.argv[2]) # size of data gap (days)
rem_width = float(sys.argv[3]) # width of data to remove from start/end of gap (days) 

gaps = []
rem = []

for i in range(1,len(time)):
	if (time[i]-time[i-1])>=gap_size:
		gaps.append(time[i-1])
		gaps.append(time[i])	

for i in range(1,len(time)):
	for j in range(len(gaps)):
		if abs(time[i]-gaps[j])<=rem_width:
			rem.append(i)

rem = list(set(rem))
rem.reverse()
new_time = np.delete(time,rem)
new_flux = np.delete(flux,rem)
new_err = np.delete(err,rem)

new_data = np.column_stack((new_time,new_flux,new_err))

np.savetxt("no_gaps.dat",new_data)
