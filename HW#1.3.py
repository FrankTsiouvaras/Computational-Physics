
#e(r) = (1/(2*math.pi^2)) integral dk k^2 P(k) sin(kr)/(kr)


# first column is k

# second column is P(k)

# r = [50,120] 

# cubic spline use mabye? .... for interpolation, cubic spline solves problem of if the function is not continuous, you cannot take the derivatives


#cubic spline vs. change of variables: (actually need both!!!!!!!!) 
##Change of variables makes the spacing uniform so Simpson's routine works,
##but it keeps the same amount of points, so the resolution is unchanged. The spline let's us
##evaluate P(k) at any k, so can have as many points as you need to resolve sin(kr).

import math
import numpy as np
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt



values = "lcdm_z0.matter_pk";

data = np.loadtxt(values);
k  = data[:,0];
Pk = data[:,1];


lnP = CubicSpline(np.log(k), np.log(Pk));

def P(kk):
    return np.exp(lnP(np.log(kk)));



# ---- Simpson on a uniform grid ----------------------------
def simpson_u(u, y):
    h = u[1] - u[0]
    return h/3*(y[0] + y[-1] + 4*y[1:-1:2].sum() + 2*y[2:-1:2].sum())

# ---- xi(r) ------------------------------------------------
def xi(r, kmin=1e-4, kmax=10.0, Npts=20001):
    u = np.linspace(np.log(kmin), np.log(kmax), Npts)   # Npts ODD -> even # intervals
    K = np.exp(u)
    y = K**3 * P(K) * np.sinc(K*r/np.pi)                # k^3 from dk = k du
    return simpson_u(u, y)/(2*np.pi**2)

# ---- evaluate and plot ------------------------------------
r  = np.linspace(50, 120, 141)
xr = np.array([xi(ri) for ri in r])

peak = r[np.argmax(r**2*xr)]

plt.plot(r, r**2*xr)
plt.axvline(peak, ls='--', c='k', lw=1)
plt.xlabel(r'$r$  [Mpc/h]')
plt.ylabel(r'$r^2\xi(r)$')
plt.title(f'BAO peak at r = {peak:.1f} Mpc/h')
plt.show()












