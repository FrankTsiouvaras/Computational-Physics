import numpy as np                      

G = 1.0                                 
M = 1.0                                 
A = 1.0                                 
B = 1.0                                 

def f(s):                               
    x = s[0]                            # pull out x
    y = s[1]                            # pull out y
    vx = s[2]                           # pull out vx
    vy = s[3]                           # pull out vy
    r = np.sqrt(x**2 + y**2)            
    v = np.sqrt(vx**2 + vy**2)          
    grav = -G*M/(4*r**3)                
    fric = -A/(v**3 + B)                
    ax = grav*x + fric*vx               
    ay = grav*y + fric*vy               
    return np.array([vx, vy, ax, ay])   

def rk4_step(s, h):                     
    k1 = h*f(s)                         
    k2 = h*f(s + 0.5*k1)                
    k3 = h*f(s + 0.5*k2)                
    k4 = h*f(s + k3)                    
    return s + (k1 + 2*k2 + 2*k3 + k4)/6  

def integrate(s, h, delta, t_end, r_stop):  
    t = 0.0                             
    tlist = [t]                         
    slist = [s]                         
    while t < t_end:                    
        r1 = rk4_step(rk4_step(s, h), h)  # way 1: two steps of size h
        r2 = rk4_step(s, 2*h)           # way 2: one step of size 2h
        diff = np.sqrt((r1[0]-r2[0])**2 + (r1[1]-r2[1])**2)  
        rho = 30*h*delta/diff           
        if rho > 1:                     
            t = t + 2*h                 
            s = r1                      
            tlist.append(t)             
            slist.append(s)             
            if np.sqrt(s[0]**2 + s[1]**2) < r_stop:  
                break                   
        h = h*min(rho**0.25, 2.0)       
    return np.array(tlist), np.array(slist)  

import matplotlib.pyplot as plt         

# ------------------- Part (a) ---------------------------
A = 0.0                                 
rp = 1e-7                               
GMp = G*M/4                             
v0 = np.sqrt(2*GMp*rp/(1 + rp))         
a_orb = (1 + rp)/2                      
T = 2*np.pi*np.sqrt(a_orb**3/GMp)       
print("v0 =", v0, " period T =", T)     

for delta in [1e-6, 1e-7, 1e-8, 1e-9, 1e-10, 1e-11]:   
    s0 = np.array([1.0, 0.0, 0.0, v0])  
    t, s = integrate(s0, 1e-3, delta, 10*T, 0.0)        # run for 10 orbits, no stopping radius
    r = np.sqrt(s[:, 0]**2 + s[:, 1]**2)                # distance at every step
    E = 0.5*(s[:, 2]**2 + s[:, 3]**2) - GMp/r           
    dE = np.max(np.abs(E - E[0]))/abs(E[0])             
    print(delta, len(t), dE, r.min())   # tolerance, number of steps, energy error, closest approach

delta = 1e-10                           
s0 = np.array([1.0, 0.0, 0.0, v0])      
t, s = integrate(s0, 1e-3, delta, 10*T, 0.0)            # run 10 orbits with the chosen tolerance
r = np.sqrt(s[:, 0]**2 + s[:, 1]**2)    # distance at every step
E = 0.5*(s[:, 2]**2 + s[:, 3]**2) - GMp/r               

plt.figure()                            
plt.plot(t/T, np.log10(r))              
plt.axhline(np.log10(rp), color="k", linestyle="--")    
plt.xlabel("t / T (orbits)")            
plt.ylabel("log10 r")                   
plt.minorticks_on()                     
plt.savefig("parta_logr.png", dpi=150)  

plt.figure()                            
plt.plot(t/T, (E - E[0])/abs(E[0]))     
plt.xlabel("t / T (orbits)")            
plt.ylabel("(E - E0) / |E0|")           
plt.minorticks_on()                     
plt.savefig("parta_energy.png", dpi=150)  

# ------------- Part (b) -----------------------------------------
A = 1.0                                 
B = 1.0                                 
vc = np.sqrt(G*M/(4*1.0))               
s0 = np.array([1.0, 0.0, 0.0, 0.8*vc])  
t, s = integrate(s0, 1e-3, delta, 1e4, 1e-7)            
r = np.sqrt(s[:, 0]**2 + s[:, 1]**2)    
print("part b: reached r =", r[-1], "at t =", t[-1], "after", len(t), "steps")  # when it got there

plt.figure()                            
plt.plot(s[:, 0], s[:, 1])              
plt.xlabel("x (100 pc)")                
plt.ylabel("y (100 pc)")                
plt.axis("equal")                       
plt.minorticks_on()                     
plt.savefig("partb_path.png", dpi=150)  

plt.figure()                            
plt.plot(t, np.log10(r))                
plt.xlabel("t (code units)")            
plt.ylabel("log10 r")                   
plt.minorticks_on()                     
plt.savefig("partb_logr.png", dpi=150)  

# ------------------- Part (c) ---------------------
# Unit conversion: r_s = 2GM/c^2 = 1e-7 x 100 pc gives M = 1.04e8 Msun,
# so one code time unit = sqrt((100 pc)^3 / (G M)) = 1.458 Myr
t_unit_Myr = 1.458                      # Myr per code time unit
delta_c = 1e-8                          # looser tolerance for the sweep (checked: merger time changes by < 2e-5 vs 1e-9, 1e-10)

pairs = [(10, 0.5), (10, 1), (5, 1), (10, 2), (4, 1), (2, 0.6), (2, 1), (10, 5), (1, 0.7),
         (0.5, 0.5), (1, 1), (5, 5), (10, 10), (1, 1.5), (0.5, 1), (1, 2), (5, 10),
         (1, 3), (1, 5), (1, 10), (0.5, 10)]  # (A, B) pairs covering B/A from 0.05 to 20
ratio = []                              
tmerge = []                             
for pair in pairs:                      # loop over the (A, B) pairs (the B/A = 10 and 20 runs take several minutes each)
    A = pair[0]                         
    B = pair[1]                         
    s0 = np.array([1.0, 0.0, 0.0, 0.8*vc])               # same start as part (b)
    t, s = integrate(s0, 1e-3, delta_c, 1e5, 1e-7)       # run until r < 1e-7
    ratio.append(B/A)                   
    tmerge.append(t[-1]*t_unit_Myr)     
    print("A =", A, " B =", B, " B/A =", B/A, " t =", t[-1]*t_unit_Myr, "Myr")  # show each result

plt.figure()                            
plt.plot(ratio, tmerge, "o")            
plt.xscale("log")                       
plt.yscale("log")                       
plt.xlabel("B / A")                     
plt.ylabel("time to reach r_s (Myr)")   
plt.savefig("partc.png", dpi=150)       

# ------------------- Part (d) ----------------------
for pair in [(1, 1), (1, 2)]:           
    A = pair[0]                         
    B = pair[1]                         
    for frac in [0.6, 0.7, 0.8, 0.9, 1.0]:                # starting speed as a fraction of circular speed
        s0 = np.array([1.0, 0.0, 0.0, frac*vc])          # start with this speed
        t, s = integrate(s0, 1e-3, delta_c, 1e5, 1e-7)   # run until r < 1e-7
        print("A =", A, " B =", B, " v0/vc =", frac, " t =", t[-1]*t_unit_Myr, "Myr")  
