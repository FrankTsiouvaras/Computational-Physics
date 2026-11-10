

# dx/dt = 1 - (b + 1)x +ax^2y ; dy/dt = bx-ax^2y
#x & y are chemical concentrations
# a & b > 0
# a = 1 ;  b = 3
# initial conditions x=y=0
# delta = 10^(-10) per unit of time in both x and y
# t=0 to t=20 ; H=20
# max n=8 (modified midpoint steps)
# if over divide in half and do it again 


import numpy as np                      

a = 1.0                                 
b = 3.0                                 
delta = 1e-10                           

def f(r):                               
    x = r[0]                            # pull out x
    y = r[1]                            # pull out y
    fx = 1 - (b+1)*x + a*x**2*y         
    fy = b*x - a*x**2*y                 
    return np.array([fx, fy])           # return both slopes as one array

def bs_interval(r, H):                  # try Bulirsch-Stoer on one interval of length H, starting from r
    
    n = 1                               
    r1 = r + 0.5*H*f(r)                 # Euler half step to the midpoint (h = H when n = 1)
    r2 = r + H*f(r1)                    # leapfrog to the end, using the slope at the midpoint
    R1 = np.empty([1, 2])               
    R1[0] = 0.5*(r1 + r2 + 0.5*H*f(r2)) 

    # ----- n = 2 to 8 -----
    while n < 8:                        
        n += 1                          
        h = H/n                         

        # modified midpoint with n steps
        r1 = r + 0.5*h*f(r)             # Euler half step to the first half point
        r2 = r + h*f(r1)                # leapfrog to the first full point
        for i in range(n-1):            
            r1 += h*f(r2)               
            r2 += h*f(r1)               

        # new row of the extrapolation table
        R2 = R1                         # previous row becomes R2
        R1 = np.empty([n, 2])           # new row has n entries
        R1[0] = 0.5*(r1 + r2 + 0.5*h*f(r2))  
        for m in range(1, n):           
            epsilon = (R1[m-1] - R2[m-1])/((n/(n-1))**(2*m) - 1)  
            R1[m] = R1[m-1] + epsilon   

        # check accuracy
        error = max(abs(epsilon[0]), abs(epsilon[1]))  
        if error < H*delta:             
            return True, R1[n-1]        

    return False, None                  

# ----- recursion -----
tpoints = [0.0]                         
xpoints = [0.0]                         
ypoints = [0.0]                         

def solve(r, t, H):                     # solve from time t to t+H, starting from r; returns r at t+H
    ok, r_new = bs_interval(r, H)       # try the whole interval in one go
    if ok:                              
        tpoints.append(t + H)           
        xpoints.append(r_new[0])        
        ypoints.append(r_new[1])        
        return r_new                    
    else:                               # it failed: split in half
        r_mid = solve(r, t, H/2)        
        return solve(r_mid, t + H/2, H/2)  # solve the second half starting from the midpoint

r_end = solve(np.array([0.0, 0.0]), 0.0, 20.0)  

# ----- plot -----
import matplotlib.pyplot as plt                       
plt.plot(tpoints, xpoints, "b-", label="x")           
plt.plot(tpoints, ypoints, "r-", label="y")           
plt.plot(tpoints, xpoints, "b.")                      
plt.plot(tpoints, ypoints, "r.")                      
plt.xlabel("t")                                       
plt.ylabel("concentration")                           
plt.legend()                                          
plt.minorticks_on()                                   
plt.savefig("brusselator.png", dpi=150)               
plt.show()                                            
