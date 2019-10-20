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

# Recalculate SNR using BLS
out_flux = v[np.argwhere(abs(phase) >= dur)]
in_flux = v[np.argwhere(abs(phase) <= 0.5*dur)]
MAD = np.median(abs(out_flux - np.median(out_flux)))
std = 1.48*MAD
best_dep = -np.mean(in_flux)
n = len(in_flux)
SNR = best_dep/std*np.sqrt(n)

# Mark all transit times
count = t0
times = []
while (count < time[-1]):
    times.append(count)
    count += per

# Test individual transits
N = 0.0
SNRm = 0.0
count = 0.0
gap = 0.0
flags = []

for i in range(len(times)):
    transit = v[np.argwhere(abs(time - times[i]) <= 0.5*dur)].flatten()
    num = len(transit)
    if num == 0:
        gap += 1.
        count += 1.
    else:
        N += 1.
        depi = -np.mean(transit)
        SNRi = depi/std*np.sqrt(num)
        if SNRi > SNRm:
            SNRm = SNRi
        exterior = time[np.argwhere(abs(time - times[i]) <= dur)]
        num_cad = len(exterior)
        num_cad_exp = 2.*dur*24.*60./29.42
        rubble = num_cad/num_cad_exp
        if (rubble < 0.75):
            gap += 1.
            flag.append(i)
        if (SNRi < 0):
            flag.append(i)

if N <= 5.:
    ses = np.zeros(len(time))
    for i in range(len(time)):
        bin_flux = v[np.argwhere(abs(time - time[i]) <= 0.5*dur)]
        ses[i] = -np.sqrt(len(bin_flux))*np.mean(bin_flux)/std
    C = []
    for i in range(len(times)):
        transit_ses = ses[np.argwhere(abs(time - times[i]) <= 0.5*dur)]
        if (len(transit_ses) > 0):
            SES_max = max(transit_ses)[0]
            features = time[np.argwhere((abs(ses) > 0.6*SES_max) & (abs(time - times[i]) >= 1.5*dur) & (abs(time - times[i]) <= per/10.))]
            if (len(features) == 0):
                C.append(1.0)
            else:
                Dt = min(abs(features - times[i]))[0]
                C.append((Dt-1.5*dur)/(per/10.))
                if (Dt-1.5*dur)/(per/10.) < 0.01:
                    flag.append(i)


list6 = list(set(flag))
count6 = count+len(list6)
N6 = len(times) - count6

gapped = gap/len(times)

# Re-calculate SNR
if len(list6) > 0:
    indices = np.argwhere(abs(time - times[list6[0]]) <= 0.5*dur)
    for i in range(1,len(list6)):
        indices = np.append(indices,np.argwhere(abs(time - times[list6[i]]) <= 0.5*dur))
    if len(indices) > 0:
        indices = sorted(indices,reverse=True)
        time2 = np.delete(time,indices)
        flux2 = np.delete(v,indices)
        phase = (((time2-t0)/per-0.5) % 1) - 0.5
        phase = phase*per
        out_transit = flux2[np.argwhere(abs(phase) >= dur)]
        in_transit = flux2[np.argwhere(abs(phase) <= 0.5*dur)]
        dep = np.mean(in_transit)
        num = len(in_transit)
        MAD = np.median(abs(out_transit - np.median(out_transit)))
        std = 1.48*MAD
        SNR6 = -np.sqrt(num)*dep/std
    else:
        SNR6 = SNR
else:
    SNR6 = SNR

if N <= 5.:
    print best_dep, SNRm, SNRm/SNR, gapped, int(N6), np.median(C), SNR6
else:
    print best_dep, SNRm, SNRm/SNR, gapped, int(N6), float('nan'), SNR6

# Returns depth of transit, max individual transit SNR, ratio of max SNR to SNR, fraction of gapped events, number of transits, re-calculated number of transits after removing flagged transits, median of chases metrics, re-calculated SNR after removing flagged transits.


	
