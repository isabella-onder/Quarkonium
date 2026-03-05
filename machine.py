from quarkonium import energy_range_finder
from quarkonium import output
from quarkonium import plotter_and_normaliser
from quarkonium import sch_solver
import constants as c
import matplotlib.pyplot as plt
import scipy as sp
import numpy as np

rmax = 30 #because it has been sufficient for now

#explicitely write and try energy range dictionary - arbitrary ranges (such that they converge)
#(could make it like my hydrogen though that would assume known peaks)
energy_ranges_c = {
    (1,0): [0,0,0.5], 
    (1,1): [0,0,1],
    (1,2): [0.8, 0, 1.5],
    (2,0): [0.5, 0, 1.3]
}

energy_ranges_b = {
    (1,0): [1, 0, 1.5],
    (1,1): [1.3,0,1.6], 
    (1,2): [1.6, 0, 1.9],
    (2,0): [1.4, 0, 1.7],
    (2,1): [1.7,0,2.0],
    (3,0): [1.8, 0, 2.2],
    (3,1): [2.1, 0, 2.3],
    (4,0): [2.1, 0, 2.3]


}

#need to update these real values with proper data
real_energy_c = {
    (1,0): 0.4380, 
    (1,1): 0.9,
    (1,2): 1.2,
    (2,0): 1.1
}


def get_energy_range(n, l, type):
    # Lookup the energy range for the given n and l
    try:
        if type == 'CHARM':
            return energy_ranges_c[(n, l)]
        elif type == 'BOTTOM':
            return energy_ranges_b[(n,l)]
    except KeyError:
        raise ValueError(f"No energy range found for n={n}, l={l}. Please check the input.")
    

def charmonium(n,l):
    quark = 'CHARM'
    E_n, _ , u, v, r =  output(l,c.m_c,c.m_c,get_energy_range(n,l,quark), c.alpha_c, c.beta_c, rmax)  #returns the final E_n and corresponding final_node (i.e. where it was cut off I think)

    print(f'\n##############################################################################\n {quark} QUARKONIUM \n The final estimated energy for n = {n}, l = {l} is E_nl = ', E_n, 'GeV', '\n############################################################################## \n' )
    return E_n, u, r

def bottonium(n,l):
    quark = 'BOTTOM'

    E_n, _ , u , v, r = output(l, c.m_b, c.m_b, get_energy_range(n,l,quark), c.alpha_b, c.beta_b, rmax )
    print(f'\n##############################################################################\n {quark} QUARKONIUM \n The final estimated energy for n = {n}, l = {l}  is E_nl = ', E_n, 'GeV', '\n############################################################################## \n' )
    return E_n, u, r

    
#charmonium(1,1)
bottonium(1,0)

def normalise(u,r):
    integral = sp.integrate.simpson(u**2,r)                           
    normalised_u = u/(np.sqrt(integral))
    return normalised_u

hyperfine = True
#making the plots
def plots():
    output_r = []
    output_u = []
    fig, axs = plt.subplots(ncols = 2, nrows = 1)
    for n_l in energy_ranges_c:
        E, u, r = charmonium(n_l[0], n_l[1])


        axs[0].scatter(r, u, label = str(n_l), marker = '.', s = 2)
        axs[0].set
        print(n_l,'was plotted successfully')
        output_r.append(r)
        output_u.append(u)
        axs[0].set_aspect(6)  
        axs[0].legend()
        axs[0].set_xlabel('Separation $r$ (GeV)')
        axs[0].set_ylabel('$u_{nl}(r)$')
    
    #not interesting to do hyperfine, radial wavefunction does not look any good
    '''if hyperfine: 
        hf_1 = 0.1604870964121239  
        _ , _, _, _, _, final_node = sch_solver(0, c.m_c, c.m_c, E + hf_1 , c.alpha_c, c.beta_c, rmax)
        nodes_nb, turning_points_nb, u_1, v, r_1, final_node = sch_solver(0, c.m_c, c.m_c, E + hf_1, c.alpha_c, c.beta_c, final_node)
        integral = sp.integrate.simpson(u_1**2,r_1)                           
        normalised_u = u_1/(np.sqrt(integral))
        #normalised_check = sp.integrate.simpson(normalised_u**2, r_1)
        axs[0].scatter(r_1, normalised_u)'''
   

    for n_l in energy_ranges_b:
        E, u, r = bottonium(n_l[0], n_l[1])
        axs[1].scatter(r, u, label = str(n_l), marker = '.', s = 2)
        print(n_l,'was plotted successfully')
        output_r.append(r)
        output_u.append(u) 
        axs[1].set_aspect(6)
        axs[1].legend()
        axs[1].set_xlabel('Separation $r$ (GeV)')
        axs[1].set_ylabel('$u_{nl}(r)$')
    plt.show()
    


#Hyperfine splitting for n = 1 is 0.1604870964121239     
#Hyperfine splitting for n = 2 is 0.11566129414755096
#plots()


