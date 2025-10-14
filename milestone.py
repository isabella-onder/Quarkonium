#%%


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#important to note: using h_bar = c = 1

def system(r, Y, l, mu, E_nl, alpha, beta): #defining my system of differential equations according to all variables
    u, v = Y                                #unpacking y 
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-alpha)/(r) + beta*r))) * u  #hmm what is beta and how do I find it + removed 4/3 for Hydrogen testing
    return [du_dr, dv_dr] 

def zero_crossing(r, Y, l, mu, E, alpha, beta):
    u, v = Y
    return u     #will trigger when the returned expression is 0

def sch_solver(l,m_1,m_2, n, alpha, beta): #for now pass everything as arguments
    mu = (1/m_1 + 1/m_2) ** (-1)
    E_nl = - mu * alpha**2/(2*n**2)
    print('this is precise E_nl', E_nl)
    E_nl = -13.602 * 10 **(-9)
    initial_conditions = [0,1] #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    r0 = 1e-6 * a0     # small start
    rmax = 10 * a0      # extend beyond peak

    print ('this is a_o', a0)
    t_eval = np.linspace(r0,rmax,1510)
    print(t_eval)

    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = t_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing) #need to find a way to put r in GeV as in get all the units right
    u, v = sol.y[0], sol.y[1]
    r = sol.t
    print('where the max occurs',r[np.argmax(u)])
    print('these are the zero crossings', sol.t_events, sol.y_events)
    print('number of crossings', len(sol.t_events[0]))

    turning_points = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] #finds turning points by checking when slope changes sign (consecutive slopes being of different sign when multiplied will be negative)
    print(r[turning_points])

    plt.scatter(r,u, marker = '.')
    plt.show()


#def energy_finder(l, m_1, m_2, E_nl, alpha, beta):
    
print(sch_solver(0,0.000511,100000000000,1,1/137,0))
    