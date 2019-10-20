import sys
import numpy as np

data = np.loadtxt(sys.argv[1])
time = data[:,0]
flux = data[:,1]
err = data[:,2]

per = float(sys.argv[2])
epo = float(sys.argv[3])
dur = float(sys.argv[4])

t0 = 54900 + epo - 0.5

count = t0
times = []

while count > time[0]:
    count -= per

while count < time[0]:
	count += per

while count > time[0] and count < time[-1]:
	times.append(count)
	count += per

a = np.where((time <= times[0]+2.0*dur) & (time >= times[0]-2.0*dur))
for i in range(1,len(times)):
	a = np.append(a,np.where((time <= times[i]+2.0*dur) & (time >= times[i]-2.0*dur)))

for i in sorted(a, reverse=True):
	time = np.delete(time, i)
	flux = np.delete(flux, i)
	err =  np.delete(err, i)


new_data = np.column_stack((time,flux,err))

np.savetxt("transit_removed.dat",new_data)
