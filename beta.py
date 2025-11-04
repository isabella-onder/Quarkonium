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


def system(r, Y, l, mu, E_nl, alpha, beta): #defining my system of differential equations: taking all parameters as input
    u, v = Y                                #unpacking y since two equations (one second order split into two first order)
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-4*alpha)/(3*r) + beta*r))) * u  #hmm what is beta and how do I find it + removed 4/3 for Hydrogen testing
    return [du_dr, dv_dr] 

def zero_crossing(r, Y, l, mu, E, alpha, beta):
    u, v = Y
    return u     #trigger function to find when u == 0, called by solve_ivp

def sch_solver(l,m_1,m_2, E_nl, alpha, beta,rmax): #passing all system parameters as arguments to make adaptable code for different particles
    mu = (1/m_1 + 1/m_2) ** (-1)
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 268101.76125244756 # Bohr radius in GeV^-1
    #a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1 that was the formula for H
    r0 = 1e-9 * a0     # small start

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


    #plt.scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
    #plt.show()

    return nodes_nb, turning_points_nb, u, r, final_node



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
            n_nb, tp_nb, _, _, final_node= sch_solver(l, m_1, m_2, energy, alpha, beta, rmax)
            nodes_nb.append(n_nb)
            turning_points_nb.append(tp_nb)
        print (nodes_nb, turning_points_nb)


        #if all nodes/turning_points are the same, means they are both on the same side as the solution: add small step to the far side and restart iteration
        #if nodes_nb[0] == nodes_nb[1] == nodes_nb[2] and turning_points_nb[0] == turning_points_nb[1] == turning_points_nb[2]: 
        #    energies[2] = energies[2] + (0.1 * 1e-9)
            #print('the energies were hopefully updated', energies)
        #    continue
        
        #replacing one side such as to narrow down depending on which side has different number of turning points/nodes
        if nodes_nb[0] != nodes_nb[1] and turning_points_nb[0] != turning_points_nb[1]:
            #print('E_1 and E_2 are different')
            betas[2] = betas[1]

        elif nodes_nb[1] != nodes_nb[2] and turning_points_nb[1] != turning_points_nb[2]:
            #print('E_2 and E_3 are different')
            betas[0] = betas[1]

        else:
                print('they are in the correct range but not yet close enough') #make sure that this really is the only scenario
                betas[0] = betas[0] + 0.0001 
                betas[2] = betas[2] - 0.0001 
        #i.e. we are not yet in the correct energy range: need to bump upwards (since that is the way we are iterating for now)
        #else:
            #print('the else loop was undergone - not yet in the correct energy range')
            #energies[0] = energies[2]
            #energies[2] = energies[2] + 0.1 * 1e-9 #the more certain we are about our range, the smaller we can make this (and if uncertain, make it larger for it to converge faster but beware if it misses it)
            #print('these are the energies after the else loop', energies)
            #continue

        

    return betas[1], final_node
    
#print(sch_solver(0,0.000511,100000000000,1,1/137,0))
#print(energy_finder(0,0.000511,100000000000,[-13.590 * 1e-9,  -13.730 * 1e-9]  ,1/137,0))


#this is just the starting point and final plotter: extracts the optimised E_nl with its final_node, and hence reruns through schrodinger equation with those parameters
#then integrates and plots
def plotter_and_normaliser(l, m_1, m_2, E_initial, alpha, beta, rmax):
    
    a0 = 268101.76125244756

    returned_betas, _ = energy_range_finder(l,m_1,m_2, E_initial, alpha, beta, rmax)

    beta = returned_betas[0] #for now just try with E_1 and hence final_node_1
    #print(final_nodes)
    #final_node = final_nodes[0]  
    rmax = 12
    print('This is rmax', rmax)

    print('This is the estimated beta', beta)
    print('these are the values going in',l, m_1, m_2, E_initial, alpha, beta, rmax )
    nodes_nb, turning_points_nb, u, r, final_node = sch_solver(l, m_1, m_2, E_initial, alpha, beta, rmax)
    print('this is nodes_nb and turning points nb of our plotted one', nodes_nb, turning_points_nb)
    #plt.scatter(r, u, marker = '.')                                  #plot the output solutions as is
    #plt.show()
  

    integral = sp.integrate.simpson(u**2,r)                           #evaluating integral over all u_nl to then normalise by result
    print('this is integral result', integral)
    normalised_u = u/(np.sqrt(integral))
    normalised_check = sp.integrate.simpson(normalised_u**2, r)
    print('this is normalised check: hopefully one', normalised_check)
    normalised_u_squared = normalised_u**2
    fig, axs = plt.subplots(1, 2)
    axs[0].scatter(r, normalised_u, marker = '.')                        #plot u_nl(r) normalised
    axs[1].scatter(r, normalised_u_squared, marker = '.')                     #plot |u_nl(r)|**2 normalised (probability density function)
    plt.show()

#plotter_and_normaliser(1,0.000511,100000000000,[-13.7 * 1e-9, 0, -13.5 * 1e-9]  ,1/137,0)
#plotter_and_normaliser(0,0.000511,100000000000,[-0.3 * 1e-9, 0, -0.2 * 1e-9]  ,1/137,0)


#attempt to get it to loop
#want to make it such that it finds all the energy levels on its own
#hence will give random lower bound and once it calculates first energy, will go upwards from there

def energy_range_finder(l,m_1,m_2, E, alpha, beta_range, rmax):                                                                              
    #energy_range = E_range
    returned_energies = []
    final_nodes = []
    for n in range(1, 2):
        beta, final_node = energy_finder(l, m_1, m_2, E, alpha, beta_range, rmax)
        print(beta, 'hopefully the final beta')
        returned_energies.append(beta)
        final_nodes.append(final_node)
        print('this is iteration n', n)
        #energy_range = [E_n + (0.1 * 1e-9), 0, E_n + (0.1 * 1e-9)]
        print('and this is the next energy range to try', beta)
    return returned_energies, final_nodes

        
plotter_and_normaliser(0,1.34,1.34, 0.388 ,0.40,[0.10,0,0.20], 10)

#M = 3.068

#energy_range_finder(0,0.000511,100000000000, 1, -13.7*1e-9,1/137,0)
#def plotter_and_normaliser(l, m_1, m_2, E_range, alpha, beta, rmax):
#plotter_and_normaliser(1,1.34 ,1.34 , , 0.40, 0, 3, 40215264.187867135)
#40215264.187867135

#need to take all the bits of final  milestone and just make sure that for quarkonium, it does not 
#require previous knowledhe om number of turns and noes (so keen the energy cjoice as anove nit omplement the rest)
