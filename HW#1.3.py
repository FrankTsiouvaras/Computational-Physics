

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


data = np.loadtxt('powerspectrum.dat')
k  = data[:,0]
Pk = data[:,1]

