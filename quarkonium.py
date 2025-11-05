#%%


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#important to note: using h_bar = c = 1



#unit converters:
def kg_to_GeVc2(m):                 #input mass in kg
    return m * 5.610 * 10 **26


def m_to_GeV(r):                    #input length in m
    return r * 5.0677 * 10 ** 15




def hydrogen_energies(n):               #for convenience: exact energy calculator for Hydrogen given n
    mu = 0.000511 #electron mass in GeV/c^-2
    alpha = 1/137 
    E_n = - mu * alpha**2/(2*n**2) #E_n in GeV
    return E_n


#defining my system of differential equations: taking all parameters as input
#unpacking y since two equations (one second order split into two first order)
def system(r, Y, l, mu, E_nl, alpha, beta): 
    u, v = Y                                
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-4*alpha)/(3*r) + beta*r))) * u  
    return [du_dr, dv_dr] 

#trigger function to find when u == 0, called by solve_ivp
def zero_crossing(r, Y, l, mu, E, alpha, beta):
    u, v = Y
    return u     

#passing all system parameters as arguments to make adaptable code for different particles
def sch_solver(l,m_1,m_2, E_nl, alpha, beta,rmax): 
    mu = (1/m_1 + 1/m_2) ** (-1)
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 268101.76125244756 # Bohr radius in GeV^-1
    r0 = 1e-10 * a0     # small start

    r_eval = np.linspace(r0,rmax,1510)  #points to evaluate u(r) at, called by solve_ivp

    #scipy function to solve differential equations system. Unpack solutions both for u and v, and corresponding distances evaluated at
    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing) 
    u, v = sol.y[0], sol.y[1] 
    r = sol.t

    #automatically finds steps where evaluates to 0 (thanks to trigger aboce)
    nodes_location = sol.t_events[0]
    print('these are the nodes_location', nodes_location)
    nodes_nb = len(nodes_location)

    #extract the final node: starts diverging from there --- NEED TO CHECK THIS QUITE BADLY
    final_node = nodes_location[-1] 
    

    #finds turning points by looking for slope change, by calculating for all pairs of points and looking for sign change when multiplied
    turning_points_index = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] 
    turning_points_location = r[turning_points_index]
    print('these are the turning points location', turning_points_location)
    turning_points_nb = len(turning_points_location)



    return nodes_nb, turning_points_nb, u, r, final_node



def energy_finder(l, m_1, m_2, energies, alpha, beta, rmax):   #input list with energy range boundaries within which to search
    #energies[2] = energies[2] + (0.1 * 1e-9)
    #print('hopefully, it has either begun or undergone a break')
    while abs(energies[0]-energies[2]) > 0.0001:

        E_2 = (energies[0] + energies[-1])/2            #bisecting energy range to start iterating
        energies = [energies[0], E_2 , energies[-1]]
        print(energies)
    

        turning_points_nb = []
        nodes_nb = []  

        #loop over 1,2,3 and store nb of nodes and turning points for each energy, in arrays where index 0 is E_1, 1 is E_2, 2 is E_3
        for E_nl in energies:  
            n_nb, tp_nb, _, _, final_node= sch_solver(l, m_1, m_2, E_nl, alpha, beta, rmax)
            nodes_nb.append(n_nb)
            turning_points_nb.append(tp_nb)
        print (nodes_nb, turning_points_nb)

        
        #replacing one side such as to narrow down depending on which side has different number of turning points/nodes
        if nodes_nb[0] != nodes_nb[1] and turning_points_nb[0] != turning_points_nb[1]:

            energies[2] = energies[1]
        elif nodes_nb[1] != nodes_nb[2] and turning_points_nb[1] != turning_points_nb[2]:
            #print('E_2 and E_3 are different')
            energies[0] = energies[1]

        else:
                print('something is wrong: the range is incorrect (or perhaps they are not yet close enough?)') #make sure that this really is the only scenario
                return('INVALID E_N', 'SO NO NODE')
                break

    

        

    return energies[1], final_node
    



#this is just the starting point and final plotter: extracts the optimised E_nl with its final_node, and hence reruns through schrodinger equation with those parameters
#then integrates and plots
def plotter_and_normaliser(l, m_1, m_2, E_initial, alpha, beta, rmax):
    a0 = 268101.76125244756
    E_nl, rmax = energy_range_finder(l,m_1,m_2, E_initial, alpha, beta, rmax)
    print('This is rmax', rmax)
    print('This is the estimated E_nl', E_nl)

    #plotting it one final time with the estimated E_nl and with rmax
    nodes_nb, turning_points_nb, u, r, final_node = sch_solver(l, m_1, m_2, E_nl, alpha, beta, rmax)
    print('this is nodes_nb and turning points nb of our plotted one', nodes_nb, turning_points_nb)

  
    #evaluating integral over all u_nl to then normalise by result
    integral = sp.integrate.simpson(u**2,r)                           
    print('this is integral result', integral)
    normalised_u = u/(np.sqrt(integral))
    normalised_check = sp.integrate.simpson(normalised_u**2, r)
    print('this is normalised check: hopefully one', normalised_check)
    normalised_u_squared = normalised_u**2

    fig, axs = plt.subplots(1, 2)
    axs[0].scatter(r/a0, normalised_u, marker = '.')                             #plot u_nl(r) normalised
    axs[1].scatter(r/a0, normalised_u_squared, marker = '.')                     #plot |u_nl(r)|**2 normalised (probability density function)
    plt.show()



#just the function that does all the calling and final returning
def energy_range_finder(l,m_1,m_2, E_range, alpha, beta, rmax):                                                                              
    energy_range = E_range
    E_n, final_node = energy_finder(l, m_1, m_2, energy_range, alpha, beta, rmax)
    print(E_n, 'hopefully the final energy')    
    print('and this is the next energy range to try', energy_range)
    return E_n, final_node

        
plotter_and_normaliser(0,1.34,1.34,[-1, 0,1],0.40,0.195,30)

