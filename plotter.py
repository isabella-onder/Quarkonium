import scipy as sp
import numpy as np
import matplotlib.pyplot as plt

#important to note: using h_bar = c = 1


def hydrogen_energies(n):               #for convenience: exact energy calculator for Hydrogen given n
    mu = 0.000511 #electron mass in GeV/c^-2
    alpha = 1/137 
    E_n = - mu * alpha**2/(2*n**2) #E_n in GeV
    print('the requested E_n', E_n)
    return E_n
print(hydrogen_energies(1), hydrogen_energies(2), hydrogen_energies(3), hydrogen_energies(4), hydrogen_energies(5))



def system(r, Y, l, mu, E_nl, alpha, beta): #defining my system of differential equations: taking all parameters as input
    u, v = Y                                #unpacking y since two equations (one second order split into two first order)
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-alpha)/(r) + beta*r))) * u  #hmm what is beta and how do I find it + removed 4/3 for Hydrogen testing
    return [du_dr, dv_dr] 

def zero_crossing(r, Y, l, mu, E, alpha, beta):
    u, v = Y
    return u     #trigger function to find when u == 0, called by solve_ivp

def sch_solver(l,m_1,m_2, n, alpha, beta): #passing all system parameters as arguments to make adaptable code for different particles
    E_nl = hydrogen_energies(n)
    E_nl = -1.37e-08
    #mu = (1/m_1 + 1/m_2) ** (-1)
    mu = 0.0005110362180000001
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    print(a0, 'this is a0')
    r0 = 1e-6 * a0     # small start
    rmax = 70* a0     # extend beyond peak
    #rmax = 7678791.00312025
    print(rmax)

    r_eval = np.linspace(r0,rmax,1510)  #points to evaluate u(r) at, called by solve_ivp

    #scipy function to solve differential equations system. Unpack solutions both for u and v, and corresponding distances evaluated at
    sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing) 
    u, v = sol.y[0], sol.y[1] 
    r = sol.t

    #automatically finds steps where evaluates to 0 thanks to trigger above)
    nodes_location = sol.t_events[0]
    nodes_nb = len(nodes_location)
    

    #finds turning points by looking for slope change, by calculating for all pairs of points and looking for sign change when multiplied
    turning_points_index = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] 
    turning_points_location = r[turning_points_index]
    turning_points_nb = len(turning_points_location)
    print(turning_points_nb, 'this is turning_points_nb')
    print('these are the turning points location', turning_points_location/a0)

    print(nodes_nb, 'this is nodes nb')
    print('these are the nodes_location', nodes_location/a0)




    fig, axs = plt.subplots(1,2)

    #axs[0].scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
    integral = sp.integrate.simpson(u**2,r)                           #evaluating integral over all u_nl to then normalise by result
    normalised_u = u/(np.sqrt(integral))
    normalised_check = sp.integrate.simpson(normalised_u**2, r)
    print('this is normalised check: hopefully one', normalised_check)

    axs[0].scatter(r/a0, normalised_u, marker = '.')                        #plot u_nl(r) normalised
    axs[1].scatter(r/a0, normalised_u**2, marker = '.')                     #plot |u_nl(r)|**2 normalised (probability density function)
    plt.show()

    return(nodes_nb, turning_points_nb, u, r)

#sch_solver(0,0.000511,100000000000,5,1/137,0)

def plotter(l,m_1,m_2, n, alpha, beta): #passing all system parameters as arguments to make adaptable code for different particles
    E_nl = hydrogen_energies(n)
    E_up = E_nl - 0.1 * 1e-9
    E_down = E_nl + 0.1 * 1e-9
    E_list = [E_down,E_nl,E_up]
    print('this is E_list', E_list)
    mu = 0.0005110362180000001
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    print(a0, 'this is a0')
    r0 = 1e-6 * a0     # small start
    rmax = 33* a0     # extend beyond peak
    #rmax = 7678791.00312025
    print(rmax)

    r_eval = np.linspace(r0,rmax,1510)  #points to evaluate u(r) at, called by solve_ivp

    fig, ax = plt.subplots()
    labels = ['+0.1eV', '$E_{3,0}$ (analytic)', '-0.1eV'] #too high, perfect, too low
    colours = ['darkseagreen', 'grey', 'peru']
    #scipy function to solve differential equations system. Unpack solutions both for u and v, and corresponding distances evaluated at
    for e,colour,label in zip(E_list, colours,labels):
        print('this is the current e', e)
        sol = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, e, alpha, beta), events = zero_crossing) 
        u, v = sol.y[0], sol.y[1] 
        r = sol.t

        #automatically finds steps where evaluates to 0 thanks to trigger above)
        nodes_location = sol.t_events[0]
        nodes_nb = len(nodes_location)
        

        #finds turning points by looking for slope change, by calculating for all pairs of points and looking for sign change when multiplied
        turning_points_index = [i for i in range(1, len(u)-1) if (u[i]-u[i-1])*(u[i+1]-u[i]) < 0] 
        turning_points_location = r[turning_points_index]
        turning_points_nb = len(turning_points_location)
        #print(turning_points_nb, 'this is turning_points_nb')
        #print('these are the turning points location', turning_points_location/a0)

        print(nodes_nb, 'this is nodes nb')
        #print('these are the nodes_location', nodes_location/a0)




  

        #axs[0].scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
        integral = sp.integrate.simpson(u**2,r)                           #evaluating integral over all u_nl to then normalise by result
        normalised_u = u/(np.sqrt(integral))
        normalised_check = sp.integrate.simpson(normalised_u**2, r)
        print('this is normalised check: hopefully one', normalised_check)

            
        plt.scatter(r/a0, normalised_u, marker = '.', s = 5, color = colour, label=label)                     #plot |u_nl(r)|**2 normalised (probability density function)
    ax.axhline(y=0, color='grey', linestyle='--')
    plt.xlabel('Separation ($a_0$)', size = 25)
    plt.ylabel('$|u_{nl}(r)|$', size = 25)
    plt.xticks(size = 15)
    plt.yticks(size = 15)
    plt.legend(markerscale=10, fontsize=17)
    plt.show()

    return(nodes_nb, turning_points_nb, u, r)

plotter(0,0.000511,100000000000,3,1/137,0)
