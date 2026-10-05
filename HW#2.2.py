

import math 

#PROBLEM #2 

#-------------------------PART B--------------------------------------------


# 5e^(-x) + x - 5 = 0
# binary search

def f(x):
    return 5*math.exp(-x)+x-5 


def Binary_Search(x1,x2):
    
    x_m = (1/2)*(x1+x2)
    
    if f(x1)*f(x2) > 0:
        print("The algorithm needs two inputs that give values of f(x) that have different signs.")
        return None

    while abs(x1-x2) > 10**(-6): 

        if f(x_m) == 0:
            print('x:',x_m) 
            return x_m
        
        elif f(x1)*f(x_m) > 0:
            x1=x_m
            x_m = (1/2)*(x1+x2)
            f(x_m) 

        elif f(x2)*f(x_m) > 0:
            x2=x_m
            x_m = (1/2)*(x1+x2)
            f(x_m) 

    print('x:',x_m) 
    return x_m 


x0 = Binary_Search(1,10)


h = 6.62607015e-34
c = 299792458
k_B = 1.380649e-23

b = (h*c)/(k_B*x0)

print("displacement constant:",b,'m\u22c5K') 


   

#-----------------------------PART C-------------------------------------------------


lam = 502e-9

T = b/lam


print("Surface Temperature of the Sun:",T,'K')
    
    
