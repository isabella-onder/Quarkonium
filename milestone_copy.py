#%%


import scipy as sp
import numpy as np
import matplotlib.pyplot as plt


#important to note: using h_bar = c = 1


#unit converters:
def kg_to_GeVc2(m):                 #input mass in kg
    return m * 5.610 * 10 **26
print(kg_to_GeVc2(9.10938 * 10 ** (-31)))

def m_to_GeV(r):                    #input length in m
    return r * 5.0677 * 10 ** 15



def hydrogen_energies(n):   
                #for convenience: exact energy calculator for Hydrogen given n
    mu = 0.000511 #electron mass in GeV/c^-2
    mu = (1/0.00051099895 + 1/0.93827208816)**(-1)
    alpha = 1/137 
    E_n = - mu * alpha**2/(2*n**2) #E_n in GeV
    print('this is the exact E_n',E_n)
    return E_n
#print('this is hydrogen energies ', hydrogen_energies(1), hydrogen_energies(2), hydrogen_energies(3), hydrogen_energies(4))



def system(r, Y, l, mu, E_nl, alpha, beta): #defining my system of differential equations: taking all parameters as input
    u, v = Y                                #unpacking y since two equations (one second order split into two first order)
    du_dr = v                               
    dv_dr = (l*(l+1))/(r**2) * u - (2 * mu * (E_nl - ((-alpha)/(r) + beta*r))) * u  #hmm what is beta and how do I find it + removed 4/3 for Hydrogen testing
    return [du_dr, dv_dr] 

def zero_crossing(r, Y, l, mu, E, alpha, beta):
    u, v = Y
    return u     #trigger function to find when u == 0, called by solve_ivp

def sch_solver(l,m_1,m_2, E_nl, alpha, beta,n,rmax): #passing all system parameters as arguments to make adaptable code for different particles
    #E_nl = -1.3605434566508328e-08 
    mu = (1/m_1 + 1/m_2) ** (-1)
    initial_conditions = [0,1]    #because we want u(0) = 0, du(0)/dr = v(0) = 1

    a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
    #print('this is a0', a0)
    r0 = 1e-6 * a0     # small start

    r_eval = np.linspace(r0,rmax,20100)  #points to evaluate u(r) at, called by solve_ivp

    #scipy function to solve differential equations system. Unpack solutions both for u and v, and corresponding distances evaluated at
    sol_2 = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing, method = 'RK23') 
    sol_3 = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing, method = 'RK45') 
    sol_4 = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing, method = 'DOP853') 
    sol_1 = sp.integrate.solve_ivp(system, [r0,rmax], initial_conditions, t_eval = r_eval, args = (l, mu, E_nl, alpha, beta), events = zero_crossing, method = 'BDF') 

    u_1, v = sol_1.y[0], sol_1.y[1] 
    r_1 = sol_1.t

    u_2, v = sol_2.y[0], sol_2.y[1] 
    r_2 = sol_2.t

    u_3, v = sol_3.y[0], sol_3.y[1] 
    r_3 = sol_3.t

    u_4, v = sol_4.y[0], sol_4.y[1] 
    r_4 = sol_4.t




    #plt.scatter(r,u, marker = '.')             #remove plotting for now since otherwise plots it every iteration
    #plt.show()

    total_u = [u_1, u_2, u_3, u_4]
    total_r = [r_1, r_2, r_3, r_4]
    return total_u, total_r

#finding the value at the origin to compare with analytic for the report: using the outputs from our optimal result.
def origin_finder():
    m_1 = 0.00051099895
    m_2 = 0.93827208816
    alpha = 1/137
    mu = (1/m_1 + 1/m_2) ** (-1)
    a0 = 1/(mu*alpha)

    #num
    _,_,u,v,r,final_node = sch_solver(0,0.00051099895,0.93827208816,-1.3605434566508328e-08 , 1/137, 0, 1, 2497579.15422085)
    2497579.15422085
    integral = sp.integrate.simpson(u**2,r)                           
    print('this is integral result', integral)
    normalised_u = u/(np.sqrt(integral))
    normalised_check = sp.integrate.simpson(normalised_u**2, r)
    print('this is the normalisation check, hopefully 1', normalised_check)
    normalised_v = v/(np.sqrt(integral))
    y_00 = 1/np.sqrt(4*np.pi)

    #analytic: for n = 1, l = 0 with exp = 0 since r = 0
    Psi = 1/(np.sqrt(np.pi)*a0**(3/2))

    print('this is normalised v at the origin', normalised_v[0])
    print('this is the num wavefunction at the origin', normalised_v[0]*y_00)
    print('this is the analytic wavefunction at the origin', Psi)
#sensitivity-wise: we have an initial difference of 0.01 eV in our energy range list, we want the absolute difference to be less thatn 0.0001 and we add 0.01 if they are not 
#in the correct range, 0.0001 if they are
#origin_finder()




def models():
    fig, axs = plt.subplots(nrows = 1, ncols =2, width_ratios=[5,2], figsize = (6,4))
    plt.subplots_adjust(wspace=0)

    u_tot, r_tot = sch_solver(0,0.00051099895,0.93827208816,-1.3605434566508328e-08 , 1/137, 0, 1, 3107579.15422085)
    normalised_u_tot = []

    integral = sp.integrate.simpson(u_tot[2]**2,r_tot[2]) #normalised according to our best result just to have meaning u on the y axis 

    for u, r in zip(u_tot, r_tot):                          
        normalised_u = u/(np.sqrt(integral))
        normalised_u_tot.append(normalised_u)
    u_tot = normalised_u_tot
    colours = ["#56A195", "#4E9BB7", "#D35328", '#2A69A4']
    colours = ["#FFB700", "#876EBD",  "#3788C9", "#FF4000"]
    colours = ["#FFB700","#36A8BA" , "#3D6E81", 'firebrick']

    models = ['BDF', 'RK3', 'RK5', 'RK8']
    energies_models = [-1.3599793090820318e-08, -1.3592001342773446e-08, -1.3605652465820318e-08, -1.360536437988282e-08]
    analytic = -1.3605434566508328e-08



    
    axs[1].set_xticks([])
    axs[1].yaxis.tick_right()
    axs[1].yaxis.set_label_position("right")
    axs[1].set_xlim(0, 1)

    for u, r, colour, model, E in zip(u_tot,  r_tot, colours, models, energies_models):
        axs[0].scatter(r/268082.760427755, u, marker = '.', color = colour, s = 7, label = model)
        axs[1].hlines(y=-E*10**9, xmin=0.15, xmax=0.85, linewidth = 3, color = colour)

    axs[0].set_xlim(xmin = 0, xmax = max(r_tot[0])/268082.760427755)
    axs[0].set_ylim(ymin = -0.2*max(u), ymax = max(u)+0.7*max(u))
    axs[0].hlines(y=0,xmin = 0, xmax = max(r)/268082.760427755,linewidth = 1.5, color='grey', linestyle = '--')

    axs[1].hlines(y=-analytic*10**9, xmin=0, xmax=1, linewidth = 2, color = 'black', linestyle = '--')
    #ax.set_xlim([0, 15])
    axs[0].set_xlabel(r'Separation $r \ (a_0)$ ', size = 12)
    axs[0].set_ylabel(r'$u_{10}(r) \ (10^{-3})$', size = 12)

    axs[1].set_ylabel('Binding energy $E_{10}$ (eV)', size = 12, labelpad = 12)

    axs[0].tick_params( direction = 'in', width = 1.4, length = 6)
    axs[1].tick_params( direction = 'in', width = 1.4, length = 6)
    axs[0].minorticks_on()
    axs[1].minorticks_on()
    axs[0].tick_params( which =  'minor', direction = 'in') 
    axs[1].tick_params( which =  'minor', direction = 'in') 
    axs[0].tick_params(which = 'major', bottom = True, left = True,  direction = 'in') 
    axs[1].tick_params(which = 'major', bottom = True, left = False,  direction = 'in') 

    axs[0].legend(markerscale=8, fontsize=12)
    from matplotlib.ticker import FuncFormatter
    axs[0].yaxis.set_major_formatter(FuncFormatter(lambda y, pos: f'{y*1000:g}'))

    #plt.savefig("summative/model_plots.svg", bbox_inches = 'tight')
    plt.savefig("summative/model_plots.png", bbox_inches = 'tight', dpi = 150)

    plt.show()

models()

