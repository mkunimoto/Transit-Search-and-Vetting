import sys
import numpy as np

data = np.loadtxt(sys.argv[1])
per = float(sys.argv[2])
epo = float(sys.argv[3])
dur = float(sys.argv[4])

time = data[:,0]
flux = data[:,1]

# Temporal time series
v = flux - np.mean(flux)

# Phased time series
t0 = epo + 54900 - 0.5
phase = ((((time - t0)/per - 0.5) % 1.) - 0.5)*per

# Recalculate SNR
out_flux = v[np.argwhere(abs(phase) >= dur)]
in_flux = v[np.argwhere(abs(phase) <= 0.5*dur)]
MAD = np.median(abs(out_flux - np.median(out_flux)))
std = 1.48*MAD
best_dep = -np.mean(in_flux)
n = len(in_flux)
SNR = best_dep/std*np.sqrt(n)

# Calculate SNR using median instead of mean
RS = -np.median(in_flux)/std*np.sqrt(n)

# Mark all transit times
count = t0
times = []
while (count < time[-1]):
    times.append(count)
    count += per

# Get CHISQR statistic
N = 0.0
chisqr = 0.0

for i in range(len(times)):
    transit = v[np.argwhere(abs(time - times[i]) <= 0.5*dur)].flatten()
    num = len(transit)
    if num != 0:
        N += 1.
        depi = -np.mean(transit)
        SNRi = depi/std*np.sqrt(num)
        chisqr += num/std**2.*(depi - best_dep)**2.

CHI = SNR*(chisqr/N)**(-0.5)

print SNR, RS, CHI, int(N)
