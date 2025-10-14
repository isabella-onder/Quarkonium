#this is my test file
#this is me trying a commit 
print('Hello, World!')
#this is me trying to change stuff on a branch

import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#testing by doing a simple spring SHM system

def system(t, Y, w): #defining my system of differential equations according to all variables
    y1, y2 = Y                                #unpacking y 
    dy1_dt = y2                             
    dy2_dt = - w**2 * y1
    return [dy1_dt, dy2_dt] 



def diff_solver(w): #for now pass everything as arguments
    initial_conditions = [1,0] #because we want x(0) = o and dx/dt(0) = 1 

    t_eval = np.linspace(0.1,15,151)

    sol = sp.integrate.solve_ivp(system, [0.1,15], initial_conditions, t_eval = t_eval, args = (w,)) #need to find a way to put r in GeV as in get all the units right
    x, v = sol.y[0], sol.y[1]
    r = sol.t


    plt.scatter(r,x)
    plt.show()

diff_solver(2)