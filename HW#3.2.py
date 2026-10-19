
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np 

# Problem 7.4____________________________________________________________________________

# Part A

xvalues = [] 
yvalues = []
x=0 
with open('dow.txt','r') as file:
    for line in file:
        if line.strip(): 
            value = float(line.strip())
            yvalues.append(value)
            xvalues.append(x)
            x+=1

plt.plot(xvalues, yvalues, color='blue', label='Dow-Jones Industrial Average')

plt.xlabel('Trading Days')
plt.ylabel('Dow-Jones Industrial Average')
plt.title('Average Prices on the US Stock Market (late 2006 Through 2010)')
plt.grid(True)
plt.legend()


ax = plt.gca() 
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))


plt.show()


#___________________________________________________________________________________


# Part B

rfft_coefficients = np.fft.rfft(yvalues)

#check = len(rfft_coefficients)
#check1 = len(yvalues)//2+1
#print(rfft_coefficients[10])
#print(check)
#print(check1) 
#513 values 

#__________________________________________________________________________________________

# Part C

first_part = round(0.1*len(rfft_coefficients))

new_rfft_coef = [] 

for i in range(first_part):
    new_rfft_coef.append(rfft_coefficients[i])

for i in range(len(rfft_coefficients) - first_part):
    new_rfft_coef.append(0) 

#___________________________________________________________________________________________


# Part D

irfft_values = np.fft.irfft(new_rfft_coef, n = len(yvalues))


plt.plot(xvalues, yvalues, color='blue', label='Dow-Jones Industrial Average')
plt.plot(xvalues, irfft_values, color='red', label='Inverse FFT (first 10% of coefficients)')

plt.xlabel('Trading Days')
plt.ylabel('Dow-Jones Industrial Average')
plt.title('Average Prices on the US Stock Market (late 2006 Through 2010)')
plt.grid(True)
plt.legend()


ax = plt.gca() 
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))


plt.show()


#_______________________________________________________________________________________________


# Part E


first_part_2 = round(0.02*len(rfft_coefficients))
new_rfft_coef_2 = [] 

for i in range(first_part_2):
    new_rfft_coef_2.append(rfft_coefficients[i])

for i in range(len(rfft_coefficients) - first_part_2):
    new_rfft_coef_2.append(0)


irfft_values_2 = np.fft.irfft(new_rfft_coef_2, n = len(yvalues))


plt.plot(xvalues, yvalues, color='blue', label='Dow-Jones Industrial Average')
plt.plot(xvalues, irfft_values_2, color='red', label='Inverse FFT (first 2% of coefficients)')

plt.xlabel('Trading Days')
plt.ylabel('Dow-Jones Industrial Average')
plt.title('Average Prices on the US Stock Market (late 2006 Through 2010)')
plt.grid(True)
plt.legend()


ax = plt.gca() 
ax.xaxis.set_major_locator(ticker.MultipleLocator(100))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(20))
ax.yaxis.set_major_locator(ticker.MultipleLocator(500))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(100))


plt.show()










