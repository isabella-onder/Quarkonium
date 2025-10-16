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

def sch_solver(l,m_1,m_2, E_nl, alpha, beta): #for now pass everything as arguments
    mu = (1/m_1 + 1/m_2) ** (-1)
    #E_nl = - mu * alpha**2/(2*n**2)
    #print('this is precise E_nl', E_nl)
    #E_nl = -13.602 * 10 **(-9)
    initial_conditions = [0,1] #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    r0 = 1e-6 * a0     # small start
    rmax = 10 * a0      # extend beyond peak

    #print ('this is a_o', a0)
    t_eval = np.linspace(r0,rmax,1510)
    #print(t_eval)

    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = t_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing) #need to find a way to put r in GeV as in get all the units right
    u, v = sol.y[0], sol.y[1]
    r = sol.t
    nodes_location = sol.t_events[0]
    nodes_nb = len(nodes_location)

    #print('where the max occurs',r[np.argmax(u)])
    #print('these are the zero crossings', sol.t_events, sol.y_events)
    #print('number of crossings', len(sol.t_events[0]))

    turning_points_index = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] #finds turning points by checking when slope changes sign (consecutive slopes being of different sign when multiplied will be negative)
    turning_points_location = r[turning_points_index]
    turning_points_nb = len(turning_points_location)
    #print('the turning points are at' , r[turning_points_index])

    

    #plt.scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
    #plt.show()

    return(nodes_nb, turning_points_nb, u, r)




def energy_finder(l, m_1, m_2, energies, alpha, beta):   #input list with energy bounds
    
    while abs(energies[0]-energies[2]) > 1e-9 * 0.001:

        E_2 = (energies[0] + energies[-1])/2
        energies = [energies[0], E_2 , energies[-1]]
        print(energies, E_2)
    

        turning_points_nb = []
        nodes_nb = []  

        for E_nl in energies:  #loop over 1,2,3 and store in arrays where index 0 is E_1, 1 is E_2, 2 is E_3
            n_nb, tp_nb, _, _= sch_solver(l, m_1, m_2, E_nl, alpha, beta)
            nodes_nb.append(n_nb)
            turning_points_nb.append(tp_nb)
        print (nodes_nb, turning_points_nb)


        #replacing it such as to narrow down
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

def plotter(l, m_1, m_2, energies, alpha, beta):
    energies = energy_finder(l, m_1, m_2, energies, alpha, beta)
    E_nl = energies[1]
    print('This is the estimated E_nl', E_nl)
    _, _, u, r = sch_solver(l, m_1, m_2, E_nl, alpha, beta)
    plt.scatter(r, u, marker = '.')
    plt.show()


plotter(0,0.000511,100000000000,[-13.5 * 1e-9, 0, -14 * 1e-9]  ,1/137,0)

    