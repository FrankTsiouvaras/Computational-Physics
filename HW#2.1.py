
import math
import numpy as np 

#6.11

# ------------- PART B --------------------------------------------------

# x = 1 - exp(-cx) c=2
#print out number of iterations it takes to converge to a solution accurate to 10^(-6)


# x* =x+e = x +e'/f'(x*)

#e' = (x-x')/[1-1/f'(x)]

#f' = ce^(-cx) 



def Relaxation_Method(x,c):
    f = 1 - math.exp(-c*x)
    df = c*math.exp(-c*x)
    error = (x-f)/(1-(1/df)) 
    iterations = 1 
    while abs(error)>10**(-6):
        x=f
        f = 1 - math.exp(-c*x)
        df = c*math.exp(-c*x)
        error = (x-f)/(1-(1/df)) 
        iterations +=1

    print(iterations)
    print(f)

    return iterations,f 



iterations_rm, x_rm = Relaxation_Method(1,2)




#---------------------- PART C ---------------------------------------------------


def OverRelaxation_Method(x,c,w):
    
    f = 1 - math.exp(-c*x)
    df = c*math.exp(-c*x)
    delta_x = f - x 
    new_x = x + (1+w)*delta_x
    error = (x - new_x)/(1 - 1/((1+w)*df-w))
    iterations = 1
    
    while abs(error)>10**(-6):
        x = new_x
        f = 1 - math.exp(-c*x)
        df = c*math.exp(-c*x)
        delta_x = f - x
        new_x = x + (1+w)*delta_x
        error = (x-new_x)/(1-1/((1+w)*df-w))
        iterations +=1

    print(iterations)
    print(new_x)

    return iterations,new_x

iterations_orm, x_orm = OverRelaxation_Method(1,2,0.5)

#need to check different values for w and compare number of iterations of Overrelaxation Method with Relaxation Method

lowest_count = float('inf') 
lowest = None 


for i in np.arange(0, 1.6, 0.1):
    iterations_orm_1, x_orm_1 = OverRelaxation_Method(1,2,i)
    print(round(i,2),iterations_orm_1, x_orm_1)
    if iterations_orm_1<lowest_count:
        lowest_count = iterations_orm_1
        lowest = round(i,2) 


rate = iterations_rm/lowest_count 
print('w =',lowest,"makes the Overrelaxation Method converge the fastest with",lowest_count,"iterations.")
print("Comparing it with the Relaxation Method, the Overrelaxation Method using w=",lowest,", the Overrelaxation Method converges",rate,"x as fast.")
        



#--------------------------PART D-----------------------------------------------


#Yes. Using ω < 0 (underrelaxation) gives faster convergence than ordinary relaxation when f′(x*) < 0.


#The overrelaxation update is x′ = x + (1 + ω)Δx, where Δx = f(x) − x is the step ordinary relaxation would take. For −1 < ω < 0, the multiplier 1 + ω is less than 1,
#so the method takes a shorter step than ordinary relaxation.

#When f′(x*) < 0, ordinary relaxation overshoots: each step carries x past the solution to the other side, so the iterates alternate above and below x*.
#Taking a shorter step reduces that overshoot and lands closer to x*, so the method converges faster.

#From the error relation in part (a), each iteration multiplies the error by
#(1 + ω) f′(x*) − ω = f′(x*) + ω [f′(x*) − 1]

#If f′(x*) < 0, then f′(x*) − 1 < 0, so choosing ω < 0 makes ω[f′(x*) − 1] positive. This moves the factor from the negative value f′(x*) toward zero,
#reducing its magnitude and so speeding up convergence. (Choosing ω > 0 here would make the factor more negative and slow convergence, or cause divergence.)

#The factor is exactly zero when
#ω = f′(x*) / (1 − f′(x*)),
#which is negative whenever f′(x*) < 0. For example, f′(x*) = -1/2 gives ω = -1/3





