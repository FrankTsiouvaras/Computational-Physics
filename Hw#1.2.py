
import math
import numpy as np 
import matplotlib.pyplot as plt


# midpoint rule


#trapezoid rule


#simpson's rule

#single precision


# integral from 0 to 1 for exp(-t)

# log-log plot of e as function of number of bins N


# errors scales as h^p (N^(-p))

#_________________________________________________________________________________


#Newman

# def f(x) : return x**4 ;   N =10 ; a =0 ; b = 2 ; h ... ; s = .5f(a) + .5f(b) ; for k in range(1,N) : s += f(a+kh) ; 


#--------------------------------------------------------------------------------

# Midpoint Method

# h = step size 

#N = (b-a)/h 

# h = (b-a)/N 

#x_i = (i + 1/2)h ; i = 0 ... N-1 (left edge)

#m_i = a + (i-1/2)*h

#I = h sum_1^N f(m_i) 

# error_t scales as h^2 (N^(-2))

#user defined input 1: N

#user-defined input 2: error tolerance: e_tol .................. (|I_(2N) - I_N|/|I_(2N)|) < e_tol



def f(x):
    return np.float32(np.exp(-np.float32(x)));



def Midpoint(a,b,N):

    a = np.float32(a);
    b = np.float32(b);
    h = np.float32((b-a)/N);
    midpoint = np.float32(0.0); 

    for i in range(1,N+1):
        m_i = a + np.float32((i -1/2))*h; 
        midpoint+= f(m_i);

    midpoint = h*midpoint;

    return midpoint; 



#-----------------------------------------------------------------------------------------


#Trapezoid Rule


# h = (b-a)/N

#x_i = ih+a


# I(a,b) = h[1/2f(a) + 1/2f(b) + sum_1^(N-1) f(x_i)]

# second order integration method (e_t scales as h^2)


def Trapezoid(a,b,N):
     
    a = np.float32(a);
    b = np.float32(b);
    h = np.float32((b-a)/N);

    trapezoid = np.float32((1/2))*f(a) + np.float32((1/2))*f(b);

    for k in range(1,N):
        x_k = k*h + a;
        trapezoid += f(x_k);

    trapezoid = h*trapezoid;

    return trapezoid; 
    


#---------------------------------------------------------------------------------------------------


#Simpson's Rule


# I(a,b)= (1/3)h[f(a)+f(b)+4sum_1^(N-1)[f(a+ih)] + 2[sum_2^(N-2)[f(a+ih)]] ....... last two summations : step size = +2


# error: 4th order (e_total scales as h^4 (N^(-4))

def Simpson(a,b,N):
    
    if N%2 != 0:
        N+=1;
        
    a = np.float32(a);
    b = np.float32(b);
    h = np.float32((b-a)/N);

    simpson = f(a) + f(b);

    r = np.float32(0.0);
    t = np.float32(0.0); 

    for i in range(1,N,2):
        r += f(a+i*h);

    for i in range(2,N-1,2):
        t += f(a+i*h); 

    simpson += (np.float32(4)*r + np.float32(2)*t);
    simpson = np.float32((1/3))*h*simpson;

    return simpson


#-------------------------------------------------------------------------------------------------------------

# Part B 

#error |I_N - I_exact|/|I_exact|

err_midpoint = [];
err_trapezoid = [];
err_simpson = [];
Ns = []; 


I_exact = 1 - math.exp(-1);

for i in range(1,21):
    N = 2**i;
    Ns.append(N); 
    
    I_1 = float(Midpoint(0,1,N));
    I_2 = float(Trapezoid(0,1,N));
    I_3 = float(Simpson(0,1,N)); 

    err_m = abs((I_1 - I_exact)/I_exact);
    err_t = abs((I_2 - I_exact)/I_exact);
    err_s = abs((I_3 - I_exact)/I_exact);

    err_midpoint.append(err_m);
    err_trapezoid.append(err_t);
    err_simpson.append(err_s);

#Plot
    
plt.loglog(Ns, err_midpoint,  'o-', label='Midpoint');
plt.loglog(Ns, err_trapezoid, 's-', label='Trapezoid');
plt.loglog(Ns, err_simpson,   '^-', label="Simpson's");

ylim = plt.ylim();

Ns = np.array(Ns);
plt.loglog(Ns, err_midpoint[0]*(Ns/Ns[0])**-2.0, 'k--', lw=2, label=r'$N^{-2}$');
plt.loglog(Ns, err_simpson[0]*(Ns/Ns[0])**-4.0,  'k:',  lw=2, label=r'$N^{-4}$');

plt.ylim(ylim);

plt.xlabel('N'); plt.ylabel('relative error')
plt.legend();
plt.show();
























