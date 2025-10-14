import scipy as sp
import numpy as np
import matplotlib as plt

def system(r, Y, l, mu, E_nl, alpha, beta): #defining my system of differential equations according to all variables
    u, v = Y                                #unpacking y 
    du_dr = v                               
    dv_dr = l*(l+1)/r**2 * u - 2 * mu(E_nl - [-(4*alpha)/(3*r) + beta*r]) * u  #hmm what is beta and how do I find it
    return [du_dr, dv_dr] 



def sch_solver(l,m_1,m_2, E_nl, alpha, beta): #for now pass everything as arguments
    mu = (1/m_1 + 1/m_2) ** (-1)
    initial_conditions = [0,1] #because we want v(0) = 1 and u(0) = 0 

    r, Y = sp.integrate.solve_ivp(system, [0,15], initial_conditions, args = (l, mu, E_nl, alpha, beta)) #need to find a way to put r in GeV as in get all the units right
    u, v = Y

    plt.scatter(r,u)
    plt.show()



    #

    

