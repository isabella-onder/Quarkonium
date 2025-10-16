#%%


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#important to note: using h_bar = c = 1


def hydrogen_energies(n):               #for convenience: exact energy calculator for Hydrogen given n
    mu = 0.000511 #electron mass in GeV/c^-2
    alpha = 1/137 
    E_n = - mu * alpha**2/(2*n**2)
    print(E_n)


def system(r, Y, l, mu, E_nl, alpha, beta): #defining my system of differential equations: taking all parameters as input
    u, v = Y                                #unpacking y since two equations (one second order split into two first order)
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-alpha)/(r) + beta*r))) * u  #hmm what is beta and how do I find it + removed 4/3 for Hydrogen testing
    return [du_dr, dv_dr] 

def zero_crossing(r, Y, l, mu, E, alpha, beta):
    u, v = Y
    return u     #trigger function to find when u == 0, called by solve_ivp

def sch_solver(l,m_1,m_2, E_nl, alpha, beta): #passing all system parameters as arguments to make adaptable code for different particles
    mu = (1/m_1 + 1/m_2) ** (-1)
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    r0 = 1e-6 * a0     # small start
    rmax = 10 * a0     # extend beyond peak


    r_eval = np.linspace(r0,rmax,1510)  #points to evaluate u(r) at, called by solve_ivp

    #scipy function to solve differential equations system. Unpack solutions both for u and v, and corresponding distances evaluated at
    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing) 
    u, v = sol.y[0], sol.y[1] 
    r = sol.t

    #automatically finds steps where evaluates to 0 (thanks to trigger aboce)
    nodes_location = sol.t_events[0]
    nodes_nb = len(nodes_location)

    #finds turning points by looking for slope change, by calculating for all pairs of points and looking for sign change when multiplied
    turning_points_index = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] 
    turning_points_location = r[turning_points_index]
    turning_points_nb = len(turning_points_location)


    #plt.scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
    #plt.show()

    return(nodes_nb, turning_points_nb, u, r)

#sch_solver(1,0.000511,100000000000,-1.5125*10**(-9)  ,1/137,0)
#sch_solver(1,0.000511,100000000000,-1.511*10**(-9)  ,1/137,0)
#sch_solver(1,0.000511,100000000000,-1.515*10**(-9)  ,1/137,0)




def energy_finder(l, m_1, m_2, energies, alpha, beta):   #input list with energy range boundaries within which to search
    
    while abs(energies[0]-energies[2]) > 1e-9 * 0.001:


        E_2 = (energies[0] + energies[-1])/2            #bisecting energy range to start iterating
        energies = [energies[0], E_2 , energies[-1]]
        print(energies, E_2)
    

        turning_points_nb = []
        nodes_nb = []  

        for E_nl in energies:  #loop over 1,2,3 and store nb of nodes and turning points for each energy, in arrays where index 0 is E_1, 1 is E_2, 2 is E_3
            n_nb, tp_nb, _, _= sch_solver(l, m_1, m_2, E_nl, alpha, beta)
            nodes_nb.append(n_nb)
            turning_points_nb.append(tp_nb)
        print (nodes_nb, turning_points_nb)


        #replacing one side such as to narrow down depending on which side has different number of turning points/nodes
        if nodes_nb[0] != nodes_nb[1] and turning_points_nb[0] != turning_points_nb[1]:
            print('E_1 and E_2 are different')
            energies[2] = energies[1]

        elif nodes_nb[1] != nodes_nb[2] and turning_points_nb[1] != turning_points_nb[2]:
            print('E_2 and E_3 are different')
            energies[0] = energies[1]
    print(energies)
    return energies
    
#print(sch_solver(0,0.000511,100000000000,1,1/137,0))
#print(energy_finder(0,0.000511,100000000000,[-13.590 * 1e-9,  -13.730 * 1e-9]  ,1/137,0))



def plotter_and_normaliser(l, m_1, m_2, energies, alpha, beta):

    energies = energy_finder(l, m_1, m_2, energies, alpha, beta)
    E_nl = energies[1]                                                #extracting final iteration of energy_finder and taking it as E_nl to evaluate at
    print('This is the estimated E_nl', E_nl)
    _, _, u, r = sch_solver(l, m_1, m_2, E_nl, alpha, beta)

    plt.scatter(r, u, marker = '.')
    plt.show()

    integral = sp.integrate.simpson(u**2,r)                           #evaluating integral over all u_nl to then normalise by result
    print('this is integral result', integral)
    normalised_u = u/(np.sqrt(integral))
    normalised_check = sp.integrate.simpson(normalised_u**2, r)
    print('this is normalised check: hopefully one', normalised_check)

    plt.scatter(r, normalised_u, marker = '.')                        #plot u_nl(r) normalised
    plt.show()
    plt.scatter(r, normalised_u**2, marker = '.')                     #plot |u_nl(r)|**2 normalised (probability density function)
    plt.show()

plotter_and_normaliser(1,0.000511,100000000000,[-1.515 * 1e-9, 0, -1.511 * 1e-9]  ,1/137,0)
#plotter_and_normaliser(0,0.000511,100000000000,[-100 * 1e-9, 0, -13.5 * 1e-9]  ,1/137,0)

    