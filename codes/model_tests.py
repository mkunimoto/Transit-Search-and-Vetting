import sys
import numpy as np
from scipy.optimize import curve_fit

# Get data to fit
data = np.loadtxt(sys.argv[1])
time = data[:,0]
flux = data[:,1]
sig = data[:,2]
guess_per = np.float(sys.argv[2])
guess_epo = np.float(sys.argv[3])
guess_dep = np.float(sys.argv[4])
guess_dur = np.float(sys.argv[5])
g1 = np.float(sys.argv[6])
g2 = np.float(sys.argv[7])
Rstar = np.float(sys.argv[8])
logg = np.float(sys.argv[9])
Teff = np.float(sys.argv[10])
guess_p = np.sqrt(guess_dep)

# Cut out all out-of-transit data
t0 = guess_epo + 54900 - 0.5
phase = (((time - t0)/guess_per - 0.5) % 1) - 0.5
phase = phase*guess_per

out_flux = flux[np.argwhere(abs(phase) > 2.5*guess_dur)].flatten()
in_time = time[np.argwhere(abs(phase) <= 2.5*guess_dur)].flatten()
in_flux = flux[np.argwhere(abs(phase) <= 2.5*guess_dur)].flatten()
in_sig = sig[np.argwhere(abs(phase) <= 2.5*guess_dur)].flatten()
MAD = np.median(abs(out_flux - np.median(out_flux)))
std = 1.48*MAD

# Define transit model - adapted from "Ian's Astro-Python Codes" (https://people.ucsc.edu/~ianc/python/)
eps = np.finfo(float).eps
zeroval = eps*1e6

def integral_smallplanet_nonlinear(z, p, cn, lower, upper):
	lower = np.array(lower, copy=True)
	upper = np.array(upper, copy=True)
	return eval_int_at_limit(upper, cn) - eval_int_at_limit(lower, cn) 


def eval_int_at_limit(limit, cn):
	sqrtlimit = np.sqrt(limit)
	sqlimit = limit*limit
	total = 1. - cn[0] * (1. - 0.8 * sqrtlimit)
	total -= cn[1] * (1. - (2./3.) * limit)
	total -= cn[2] * (1. - (4./7.) * limit*sqrtlimit)
	total -= cn[3] * (1. - 0.5 * sqlimit)
	ret = -(sqlimit) * total
	return ret


def transit_model(time, epo, per, bp, ars, p, zpt):
    if (bp < 0.0) or (ars < 1.0) or (p < 0.0) or (bp > 1.0):
        return 1000.
    else:
        cn = [0,g1+2.*g2,0,-g2]
        z = t2z(time, epo, per, bp, ars)
        cn = np.array([cn], copy=False).ravel()
        if cn.size < 4:
            cn = np.concatenate((cn, [0.]*(4-cn.size)))
        z = np.array(z, copy=False)
        F = np.ones(z.shape, float)
        z[z==0] = zeroval
        a = (z - p)**2
        b = (z + p)**2
        c0 = 1. - np.sum(cn)
        Omega = 0.25 * c0 + np.sum( cn / np.arange(5., 9.) )
        ind1 = ((1. - p) < z) * ((1. + p) > z)
        ind2 = z <= (1. - p)
        aind1 = 1. - a[ind1]
        zind1m1 = z[ind1] - 1.
        Istar_edge = integral_smallplanet_nonlinear(None, p, cn, np.sqrt(aind1), np.array([0.])) / aind1
        Istar_inside = integral_smallplanet_nonlinear(None, p, cn, np.sqrt(1. - a[ind2]), np.sqrt(1. - b[ind2])) / (z[ind2])
        term1 = 0.25 * Istar_edge / (np.pi * Omega)
        term2 = p*p * np.arccos((zind1m1) / p)
        term3 = (zind1m1) * np.sqrt(p*p - (zind1m1*zind1m1))
        F[ind1] = 1. - term1 * (term2 - term3)
        F[ind2] = 1. - 0.0625 * p * Istar_inside / Omega
        return F-1+zpt


def t2z(time, epo, per, bp, ars):
    t0 = epo + 54900 - 0.5
    omega_tdiff = (2*np.pi/per)*(time-t0)
    cosom = np.cos(omega_tdiff)
    z = ars*np.sqrt(np.sin(omega_tdiff)**2. + ((bp/ars)*cosom)**2.)
    z[cosom<0] = 100.
    return z

# Fit transit models with different initial guess for impact parameter
bp = np.arange(0.0,1.0,0.1)

i = 0
fit = 0
model_red = np.inf
while i < len(bp):
    try:
        sterm = np.sin(np.pi*guess_dur/guess_per)
        guess_ars = ((1 + guess_p)**2. - bp[i]*bp[i]*(1 - sterm**2.))/sterm**2.
        guess_ars = np.sqrt(guess_ars)
        popt, pcov = curve_fit(transit_model, in_time, in_flux, p0 = [guess_epo, guess_per, bp[i], guess_ars, guess_p, 0.0])
        try:
            f_model = transit_model(in_time,*popt)
            chisqr = np.sum(((in_flux - f_model)/in_sig)**2.)
            red_chisqr = chisqr/(len(in_time) - 6)
            if np.abs(red_chisqr - 1) < np.abs(model_red - 1):
                fit = 1
                best_vals = popt
                model_rss = np.sum((in_flux - f_model)**2.)
                model_chi = chisqr
                model_dof = len(in_time) - 6
                model_red = red_chisqr
        except ValueError:
            cont = 0
    except RuntimeError:
        cont = 0
    i += 1

# Fit straight line model
def straight_model(t, l):
    return -l*np.ones(len(t))

popt, pcov = curve_fit(straight_model, in_time, in_flux, p0 = [0.5*guess_dep])
s = straight_model(in_time,*popt)
straight_rss = np.sum((in_flux - s)**2.)
straight_chi = np.sum(((in_flux - s)/in_sig)**2.)
straight_dof = len(in_time) - 1
straight_red = straight_chi/straight_dof

# Fit trapezoid model
def trapezoid(x, x_0, dep, width, slope, zpt):
    if (dep < 0) or (slope < 0):
        return 1000.
    else:
        x2 = x_0 - width/2.
        x3 = x_0 + width/2.
        x1 = x2 - dep/slope
        x4 = x3 + dep/slope
        range_a = np.logical_and(x >= x1, x < x2)
        range_b = np.logical_and(x >= x2, x < x3)
        range_c = np.logical_and(x >= x3, x < x4)
        range_d = np.logical_and(x >= min(x), x < x1)
        range_e = np.logical_and(x >= x4, x < max(x))
        val_a = zpt-slope*(x-x1)
        val_b = zpt-dep
        val_c = zpt-slope*(x4-x)
        val_d = zpt
        result = np.select([range_a, range_b, range_c, range_d, range_e], [val_a, val_b, val_c, val_d, val_d])
        return result

if (fit == 1):
    f_model = transit_model(in_time,*best_vals)
    best_epo = best_vals[0]
    best_per = best_vals[1]
    best_b = best_vals[2]
    best_ars = best_vals[3]
    best_p = best_vals[4]
    best_zpt = best_vals[5]
    best_dur = (best_per/np.pi)*np.arcsin(np.sqrt(((1+best_p)**2. - best_b**2.)/(best_ars**2. - best_b**2.)))
elif (fit == 0):
    t0 = guess_epo + 54900 - 0.5
    phase = (((in_time - t0)/guess_per - 0.5) % 1) - 0.5
    phase = phase*guess_per
    popt, pcov = curve_fit(trapezoid, phase, in_flux, p0 = [0.0, guess_dep, 0.8*guess_dur, guess_dep/0.1/guess_dur, 0.0])
    f_model = trapezoid(phase,*popt)
    model_rss = np.sum((in_flux - f_model)**2.)
    model_chi = np.sum(((in_flux - f_model)/in_sig)**2.)
    model_dof = len(in_time) - 5
    model_red = model_chi/model_dof
    best_epo = guess_epo + popt[0]
    best_per = guess_per
    best_dep = popt[1]
    best_p = np.sqrt(best_dep)
    best_zpt = popt[4]
    best_width = popt[2]
    best_slope = popt[3]
    best_in = best_dep/best_slope
    best_dur = best_width + 2.*best_in
    sinF = (np.sin(best_width*np.pi/best_per))**2.
    sinT = (np.sin(best_dur*np.pi/best_per))**2.
    best_b = np.sqrt(((1 - best_p)**2. - sinF/sinT*(1 + best_p)**2.)/(1 - sinF/sinT))
    best_ars = np.sqrt(((1 + best_p)**2. - best_b**2.*(1 - sinT))/sinT)


SNR = np.sqrt(np.dot(f_model-best_zpt,f_model-best_zpt))/std
vshape = best_b + best_p

Rsun = 6.95508e8 # radius of sun
Rearth = 6.371e6 # radius of earth
Tsun = 5778. # temperature of sun
fau = 1.496e11 # m per au
g = 10.**(logg)/100.

Rp = (Rsun/Rearth)*best_p*Rstar
a = ((86400*best_per)*(Rstar*Rsun)*np.sqrt(g)/(2.*np.pi))**(2./3.)/fau
Teq = Teff*(1 - 0.3)**(0.25)*np.sqrt(Rstar*Rsun/(2.*a*fau))
Seff = (Rstar/a)**2.*(Teff/Tsun)**4.

# Odd-even transits:
odd_t0 = best_epo + 54900 - 0.5
even_t0 = best_epo + 54900 - 0.5 + best_per
odd_phase = (((time - odd_t0)/(2.*best_per) - 0.5) % 1.)*2.*best_per
even_phase = (((time - even_t0)/(2.*best_per) - 0.5) % 1.)*2.*best_per

if (fit == 1):
    odd_flux = flux[np.argwhere(abs(odd_phase - best_per) <= 30./1440.)]
    even_flux = flux[np.argwhere(abs(even_phase - best_per) <= 30./1440.)]
elif (fit == 0):
    if (best_width > 1./24.):
        odd_flux = flux[np.argwhere(abs(odd_phase - best_per) <= 0.5*best_width)]
        even_flux = flux[np.argwhere(abs(even_phase - best_per) <= 0.5*best_width)]
    else:
        odd_flux = flux[np.argwhere(abs(odd_phase - best_per) <= 30./1440.)]
        even_flux = flux[np.argwhere(abs(even_phase - best_per) <= 30./1440.)]

odd_dep = np.median(odd_flux)
odd_err = np.std(odd_flux)
even_dep = np.median(even_flux)
even_err = np.std(even_flux)
oe = abs(odd_dep - even_dep)/np.sqrt(odd_err**2. + even_err**2.)


print fit, SNR, model_red, straight_red, Rp, a, Teq, Seff, best_per, best_epo, best_dur, best_p, best_b, vshape, oe

# Save best-fit model with parameters
if (fit == 1):
    f_model = transit_model(time,*best_vals)
    mat = np.column_stack((time,flux,f_model))
    np.savetxt("LMfit_model.dat",mat)
    file = open('LMfit_params.dat', 'w')
    file.write(str(1) + ' ' + str(best_per) + ' ' + str(best_epo) + ' ' + str(best_p) + ' ' + str(best_ars) + ' ' + str(best_b) + ' ' + str(best_zpt))
    file.close()
elif (fit == 0):
    f_model = trapezoid(time,*popt)
    mat = np.column_stack((time,flux,f_model))
    np.savetxt("LMfit_model.dat",mat)
    file = open('LMfit_params.dat', 'w')
    file.write(str(0) + ' ' + str(best_per) + ' ' + str(best_epo) + ' ' + str(popt[1]) + ' ' + str(best_width) + ' ' + str(best_slope) + ' ' + str(best_zpt))
    file.close()
