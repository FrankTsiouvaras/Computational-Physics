
import math
import numpy as np
from pylab import plot, show
import matplotlib.pyplot as plt

#1 Problem 1 a

# cos(x) and exp(x) at x= 0.1, 10 .... use single precision 

#----------------------------------------------------------------------------------------------
# FORWARD DIFFERENCE 

# df = [f(x+h) - f(x)]/h

# error proportional to h

# h is roughly square root of machine error, which in single precision is 10^(-7) ... np.float()

# h optimal = (e_m)^(1/2) ... 10^(-3.5)

#for e in np.arange(0,-8.2,-0.2):
#    h = 10**e;

def forward(f,x,h):
    return (f(x+h) -f(x))/h 

#-----------------------------------------------------------------------------------------

# CENTRAL DIFFERENCE

# df = [f(x+h)-f(x-h)]/(2h)

# total error scales as (machine error)^(2/3)

# h optimal = (e_m)^(1/3) ... 10^(-2.5)

def central(f,x,h):
    return (f(x+h)-f(x-h))/(np.float32(2.0)*h)

#-------------------------------------------------------------------------------------------

# EXTRAPOLATION

# 2 values of h

# h1 = h1
# h2 = 2h1

# df = [-f(x+2h) + 8f(x+h) - 8f(x-h) + f(x-2h)]/(12h)


# total error scalses as (machine error)^(4/5) !!!!!!! (notes on sept 15 say (3/5) !!!!!!!


# h optimal = (e_m)^(1/5) ... 10^(-1.5)

def extrapolation(f,x,h):
    return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/(np.float32(12)*h)

#--------------------------------------------------------------------------------------------------

#main 


x1 = np.float64(np.float32(0.1))
x2 = 10.0

exact_cos_01 = -np.sin(x1)
exact_cos_10 = -np.sin(x2)
exact_exp_01 = np.exp(x1)
exact_exp_10 = np.exp(x2)


err_cos_01_for = [];
err_cos_10_for = [];
err_exp_01_for = [];
err_exp_10_for = [];


err_cos_01_cen = [];
err_cos_10_cen = [];
err_exp_01_cen = [];
err_exp_10_cen = [];


err_cos_01_ex = [];
err_cos_10_ex = [];
err_exp_01_ex = [];
err_exp_10_ex = [];

h_values = []; 

for p in np.arange(0,-8.2,-0.2):
    h = np.float32(10**p);
    
    cos_f1 = forward(np.cos,np.float32(.1),h);
    cos_f1_e = abs((cos_f1 - exact_cos_01)/exact_cos_01); 

    
    cos_f2 = forward(np.cos,np.float32(10),h);
    cos_f2_e = abs((cos_f2 - exact_cos_10)/exact_cos_10); 

    
    ex_f1 = forward(np.exp,np.float32(.1),h);
    ex_f1_e = abs((ex_f1 - exact_exp_01)/exact_exp_01); 

    
    ex_f2 = forward(np.exp,np.float32(10),h);
    ex_f2_e = abs((ex_f2 - exact_exp_10)/exact_exp_10);

    cos_c1 = central(np.cos,np.float32(.1),h);
    cos_c1_e = abs((cos_c1 - exact_cos_01)/exact_cos_01); 

    
    cos_c2 = central(np.cos,np.float32(10),h);
    cos_c2_e = abs((cos_c2 - exact_cos_10)/exact_cos_10); 

    
    ex_c1 = central(np.exp,np.float32(.1),h);
    ex_c1_e = abs((ex_c1 - exact_exp_01)/exact_exp_01); 

    
    ex_c2 = central(np.exp,np.float32(10),h);
    ex_c2_e = abs((ex_c2 - exact_exp_10)/exact_exp_10);

    cos_expo1 = extrapolation(np.cos,np.float32(.1),h);
    cos_expo1_e = abs((cos_expo1 - exact_cos_01)/exact_cos_01); 

    
    cos_expo2 = extrapolation(np.cos,np.float32(10),h);
    cos_expo2_e = abs((cos_expo2 - exact_cos_10)/exact_cos_10); 

    
    ex_expo1 = extrapolation(np.exp,np.float32(.1),h);
    ex_expo1_e = abs((ex_expo1 - exact_exp_01)/exact_exp_01); 

    
    ex_expo2 = extrapolation(np.exp,np.float32(10),h);
    ex_expo2_e = abs((ex_expo2 - exact_exp_10)/exact_exp_10);  



    
    err_cos_01_for.append(cos_f1_e);
    err_cos_10_for.append(cos_f2_e);
    err_exp_01_for.append(ex_f1_e);
    err_exp_10_for.append(ex_f2_e);

    err_cos_01_cen.append(cos_c1_e);
    err_cos_10_cen.append(cos_c2_e);
    err_exp_01_cen.append(ex_c1_e);
    err_exp_10_cen.append(ex_c2_e);

    err_cos_01_ex.append(cos_expo1_e);
    err_cos_10_ex.append(cos_expo2_e);
    err_exp_01_ex.append(ex_expo1_e);
    err_exp_10_ex.append(ex_expo2_e);

    h_values.append(h); 



# minima


i1 = np.argmin(err_cos_01_for);
i2 = np.argmin(err_cos_01_cen);
i3 = np.argmin(err_cos_01_ex);

print("cos(x) at x = 0.1");
print("forward:", h_values[i1], err_cos_01_for[i1]);
print("central:", h_values[i2], err_cos_01_cen[i2]);
print("extrapolated:", h_values[i3], err_cos_01_ex[i3]);


#-----------------------------------------------------------------------------

i4 = np.argmin(err_cos_10_for);
i5 = np.argmin(err_cos_10_cen);
i6 = np.argmin(err_cos_10_ex);

print("cos(x) at x = 10");
print("forward:", h_values[i4], err_cos_10_for[i4]);
print("central:", h_values[i5], err_cos_10_cen[i5]);
print("extrapolated:", h_values[i6], err_cos_10_ex[i6]);


#-----------------------------------------------------------------------------

i7 = np.argmin(err_exp_01_for);
i8 = np.argmin(err_exp_01_cen);
i9 = np.argmin(err_exp_01_ex);

print("exp(x) at x = 0.1");
print("forward:", h_values[i7], err_exp_01_for[i7]);
print("central:", h_values[i8], err_exp_01_cen[i8]);
print("extrapolated:", h_values[i9], err_exp_01_ex[i9]);


#-----------------------------------------------------------------------------

i10 = np.argmin(err_exp_10_for);
i11 = np.argmin(err_exp_10_cen);
i12 = np.argmin(err_exp_10_ex);

print("exp(x) at x = 10");
print("forward:", h_values[i10], err_exp_10_for[i10]);
print("central:", h_values[i11], err_exp_10_cen[i11]);
print("extrapolated:", h_values[i12], err_exp_10_ex[i12]);

#--------------------------------------------------------------------------------

# PART B  


# relative error = |(calculated value) - (true value)|/|(true value)|


# PLOTS!!!!!!!



#-----------------------------------------------------------------------------------------------------
h_ref = np.array(h_values);
eps = np.finfo(np.float32).eps;

plt.figure();

plt.xlabel("h");
plt.ylabel("relative error");
plt.title("cos(x) at x = 0.1");

plt.loglog(h_values, err_cos_01_for, "r-", label="forward");
plt.loglog(h_values, err_cos_01_cen, "g-", label="central");
plt.loglog(h_values, err_cos_01_ex,  "c-", label="extrapolated");
plt.loglog(h_ref, eps*abs(np.cos(0.1)/np.sin(0.1))/h_ref, "k--", label=r"roundoff, $\epsilon|f/f'|/h$");
plt.legend();

plt.loglog(h_values[i1], err_cos_01_for[i1], "ro");
plt.loglog(h_values[i2], err_cos_01_cen[i2], "go");
plt.loglog(h_values[i3], err_cos_01_ex[i3],  "co");

plt.grid(False);
plt.tick_params(which="major", length=8, direction="in");
plt.tick_params(which="minor", length=4, direction="in");

#-----------------------------------------------------------------------------

plt.figure();


plt.xlabel("h");
plt.ylabel("relative error");
plt.title("cos(x) at x = 10");

plt.loglog(h_values, err_cos_10_for, "r-", label="forward");
plt.loglog(h_values, err_cos_10_cen, "g-", label="central");
plt.loglog(h_values, err_cos_10_ex,  "c-", label="extrapolated");
plt.loglog(h_ref, eps*abs(np.cos(10)/np.sin(10))/h_ref, "k--", label=r"roundoff, $\epsilon|f/f'|/h$");
plt.legend();

plt.loglog(h_values[i4], err_cos_10_for[i4], "ro");
plt.loglog(h_values[i5], err_cos_10_cen[i5], "go");
plt.loglog(h_values[i6], err_cos_10_ex[i6],  "co");

plt.grid(False);
plt.tick_params(which="major", length=8, direction="in");
plt.tick_params(which="minor", length=4, direction="in");


#---------------------------------------------------------------------------------------

plt.figure();

plt.xlabel("h");
plt.ylabel("relative error");
plt.title("exp(x) at x = 0.1");

plt.loglog(h_values, err_exp_01_for, "r-", label="forward");
plt.loglog(h_values, err_exp_01_cen, "g-", label="central");
plt.loglog(h_values, err_exp_01_ex,  "c-", label="extrapolated");
plt.loglog(h_ref, eps/h_ref, "k--", label=r"roundoff, $\epsilon|f/f'|/h$");
plt.legend();

plt.loglog(h_values[i7], err_exp_01_for[i7], "ro");
plt.loglog(h_values[i8], err_exp_01_cen[i8], "go");
plt.loglog(h_values[i9], err_exp_01_ex[i9],  "co");

plt.grid(False);
plt.tick_params(which="major", length=8, direction="in");
plt.tick_params(which="minor", length=4, direction="in");


#---------------------------------------------------------------------------------------------



plt.figure();

plt.xlabel("h");
plt.ylabel("relative error");
plt.title("exp(x) at x = 10");

plt.loglog(h_values, err_exp_10_for, "r-", label="forward");
plt.loglog(h_values, err_exp_10_cen, "g-", label="central");
plt.loglog(h_values, err_exp_10_ex,  "c-", label="extrapolated");
plt.loglog(h_ref, eps/h_ref, "k--", label=r"roundoff, $\epsilon|f/f'|/h$");
plt.legend();

plt.loglog(h_values[i10], err_exp_10_for[i10], "ro");
plt.loglog(h_values[i11], err_exp_10_cen[i11], "go");
plt.loglog(h_values[i12], err_exp_10_ex[i12],  "co");

plt.grid(False);
plt.tick_params(which="major", length=8, direction="in");
plt.tick_params(which="minor", length=4, direction="in");



plt.show();
    




