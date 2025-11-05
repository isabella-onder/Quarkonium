#%%


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#important to note: using h_bar = c = 1



#unit converters:
def kg_to_GeVc2(m):                 #input mass in kgm, returns in Gev/c^2
    return m * 5.610 * 10 **26


def m_to_GeV(r):                    #input length in m, returns in Gev^-1
    return r * 5.0677 * 10 ** 15



#defining my system of differential equations: taking all parameters as input
#unpacking y since two equations (one second order split into two first order)
def system(r, Y, l, mu, E_nl, alpha, beta): 
    u, v = Y                                
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-4*alpha)/(3*r) + beta*r))) * u 
    return [du_dr, dv_dr] 

#trigger function to find when u == 0, called by solve_ivp, to flag nodes
def zero_crossing(r, Y, l, mu, E, alpha, beta): 
    u, v = Y
    return u     


#passing all system parameters as arguments to make adaptable code for different particles
def sch_solver(l,m_1,m_2, E_nl, alpha, beta,rmax): 
    mu = (1/m_1 + 1/m_2) ** (-1)
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 268101.76125244756 # Bohr radius in GeV^-1
    r0 = 1e-9 * a0     # very small start to avoid odd behaviour at 0

    r_eval = np.linspace(r0,rmax,1510)  #points to evaluate u(r) at, called by solve_ivp

    #scipy function to solve differential equations system. Unpack solutions both for u and v, and corresponding distances evaluated at
    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing) 
    u, v = sol.y[0], sol.y[1] 
    r = sol.t

    #automatically finds steps where evaluates to 0 (thanks to trigger above)
    nodes_location = sol.t_events[0]
    print('these are the nodes_location', nodes_location)
    nodes_nb = len(nodes_location)

    #extract the final node: starts diverging from there --- DOES NOT WORK FOR QUARKONIUM IT SEEMS
    #final_node = nodes_location[-1] 
    

    #finds turning points by looking for slope change, by calculating for all pairs of points and looking for sign change when multiplied
    turning_points_index = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] 
    turning_points_location = r[turning_points_index]
    print('these are the turning points location', turning_points_location)
    turning_points_nb = len(turning_points_location)


    #plt.scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
    #plt.show()

    return nodes_nb, turning_points_nb, u, r



def energy_finder(l, m_1, m_2, energy, alpha, betas, rmax):   #input list with energy range boundaries within which to search
    #energies[2] = energies[2] + (0.1 * 1e-9)
    #print('hopefully, it has either begun or undergone a break')
    while abs(betas[0]-betas[2]) > 0.0001:

        beta_2 = (betas[0] + betas[-1])/2            #bisecting energy range to start iterating
        betas = [betas[0], beta_2 , betas[-1]]
        print(betas)
    

        turning_points_nb = []
        nodes_nb = []  

        for beta in betas:  #loop over 1,2,3 and store nb of nodes and turning points for each energy, in arrays where index 0 is E_1, 1 is E_2, 2 is E_3
            n_nb, tp_nb, _, _= sch_solver(l, m_1, m_2, energy, alpha, beta, rmax)
            nodes_nb.append(n_nb)
            turning_points_nb.append(tp_nb)
        print (nodes_nb, turning_points_nb)
        
        #replacing one side such as to narrow down depending on which side has different number of turning points/nodes
        if nodes_nb[0] != nodes_nb[1] and turning_points_nb[0] != turning_points_nb[1]:
            #print('E_1 and E_2 are different')
            betas[2] = betas[1]

        elif nodes_nb[1] != nodes_nb[2] and turning_points_nb[1] != turning_points_nb[2]:
            #print('E_2 and E_3 are different')
            betas[0] = betas[1]

        else:
                print('they are either not close enough or on the wrong side alltogether')
                return('NO VALID BETA')
                break  

        
    return betas[1]
    


#this  plotter: extracts the optimised E_nl with its final_node, and hence reruns through schrodinger equation with those parameters
#then integrates and plots
def plotter_and_normaliser(l, m_1, m_2, E_initial, alpha, beta, rmax):
    
    a0 = 268101.76125244756

    returned_betas = beta_range_finder(l,m_1,m_2, E_initial, alpha, beta, rmax)

    beta = returned_betas[0] #for now just try with E_1 and hence final_node_1
    #print(final_nodes)
    #final_node = final_nodes[0]  
    rmax = 12
    print('This is rmax', rmax)

    print('This is the estimated beta', beta)
    print('these are the values going in',l, m_1, m_2, E_initial, alpha, beta, rmax )
    nodes_nb, turning_points_nb, u, r = sch_solver(l, m_1, m_2, E_initial, alpha, beta, rmax)
    print('this is nodes_nb and turning points nb of our plotted one', nodes_nb, turning_points_nb)
 
  
    #evaluating integral over all u_nl to then normalise by result & checking it is indeed 1 afterwards
    integral = sp.integrate.simpson(u**2,r)                           
    normalised_u = u/(np.sqrt(integral))
    normalised_check = sp.integrate.simpson(normalised_u**2, r)
    print('this is normalised check: hopefully one', normalised_check)
    
    normalised_u_squared = normalised_u**2
    fig, axs = plt.subplots(1, 2)
    axs[0].scatter(r, normalised_u, marker = '.')                        #plot u_nl(r) normalised
    axs[1].scatter(r, normalised_u_squared, marker = '.')                #plot |u_nl(r)|**2 normalised (probability density function)
    plt.show()



#just the function that does all the calling and final returning
def beta_range_finder(l,m_1,m_2, E, alpha, beta_range, rmax):                                                                              
    returned_energies = []
    #final_nodes = []
    beta = energy_finder(l, m_1, m_2, E, alpha, beta_range, rmax)
    print(beta, 'hopefully the final beta')
    returned_energies.append(beta)
    return returned_energies

        
#plotter_and_normaliser(0,1.34,1.34, 0.388 ,0.40,[0.10,0,0.20], 10)
beta_range_finder(0,1.34,1.34, 0.388 ,0.40,[0.10,0,0.20], 10)

