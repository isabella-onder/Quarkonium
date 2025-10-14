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



def sch_solver(l,m_1,m_2, n, alpha, beta): #for now pass everything as arguments
    mu = (1/m_1 + 1/m_2) ** (-1)
    E_nl = - mu * alpha**2/(2*n**2)
    print(E_nl)
    #print('this is mu', mu)
    initial_conditions = [0,1] #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    r0 = 1e-6 * a0     # small start
    rmax = 20 * a0      # extend beyond peak

    print ('this is a_o', a0)
    t_eval = np.linspace(r0,rmax,1510)
    print(t_eval)

    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = t_eval, args = (l, mu, E_nl, alpha, beta)) #need to find a way to put r in GeV as in get all the units right
    u, v = sol.y[0], sol.y[1]
    r = sol.t
    print(r)

    plt.scatter(r,u)
    plt.show()


print(sch_solver(1,0.000511,100000000000,2,1/137,0))
    


# %%
