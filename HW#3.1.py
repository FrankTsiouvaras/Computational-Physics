
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker


#_________________________Problem 7.2_____________________________________________





#Part A

coordinates = [] 

with open('sunspots.txt','r') as file:
    for line in file:
        parts = line.strip()
        parts = parts.split()

        if len(parts)==2:
            x = float(parts[0])
            y = float(parts[1])

            coordinates.append((x,y))

            
x_values = []
y_values = []

for tup in coordinates:
    x_values.append(tup[0])
    y_values.append(tup[1])

    
plt.plot(x_values, y_values, color='blue', label='Sunspots')


plt.xlabel('Time in Months')
plt.ylabel('Sunspots')
plt.title('Periodicity of Sunspots through Time')
plt.grid(True)
plt.legend()


ax = plt.gca() 
ax.xaxis.set_major_locator(ticker.MultipleLocator(500))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(100))
ax.yaxis.set_major_locator(ticker.MultipleLocator(50))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(10))


plt.show()


# estimate of length of cycle in months: 130 months


#________________________________________________________________________________________

# Part B


# Fourier Coefficients : alpha_k = 1/L integral 0 to L f(x)exp(-i (2pikx/L)dx
# alpha_k = 1/N sum (0 to N-1) f(x_n)exp(-i 2pikx_n/L) 
#discrete points : x_n = (n/N)L ..... x_n/L = n/N
# alpah_k = 1/N sum (0 to N-1) y_n exp(-i 2pikn/N)
# alpha_k = C_k/N
#C_k .. 0 to N/2 
# sqare wave

import numpy as np


coordinates = [] 
N = 0

with open('sunspots.txt','r') as file:
    for line in file:
        parts = line.strip()
        parts = parts.split()

        if len(parts)==2:
            x = float(parts[0])
            y = float(parts[1])

            coordinates.append((x,y))
            N+=1

            
x_values = []
y_values = []

for tup in coordinates:
    x_values.append(tup[0])
    y_values.append(tup[1])



coefficients = []
k_values = [] 

for k in range(0,N//2+1):
    C_k = 0
    k_values.append(k) 
    for n in range(0,N):
        C_k += y_values[n]*np.exp(-2j * np.pi * k * n / N)
    coefficients.append(C_k) 

y_new_values = []
for item in coefficients:
    y_new_values.append(abs(item)**2)




plt.plot(k_values, y_new_values, color='blue', label='|C_k|^2')


plt.xlabel('k values')
plt.ylabel('|C_k|^2')
plt.title('Power Spectrum of Sunspots')
plt.grid(True)
plt.legend()

plt.ylim(0, 1.1 * max(y_new_values[1:])) # so C_0 does not dominate/squash everything
plt.xlim(0,200) #zoom in on low k area where peak is 
ax = plt.gca()
ax.xaxis.set_major_locator(ticker.MultipleLocator(20))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))

plt.show()
#slower runtime that is should probably be ..............



#_________________________________________________________________________________________________


# Part C


ymax_value = 0
kmax_value = 0 
for i in range(1, N//2+1):
    if y_new_values[i] > ymax_value:
        ymax_value = y_new_values[i] 
        kmax_value = k_values[i] 
    
print(kmax_value)

# k cycles fit into N months 
period = N/kmax_value

print(period) 

# matches pretty closley to approx in part (a) 
