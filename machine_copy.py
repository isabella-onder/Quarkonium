from quarkonium import energy_range_finder
from quarkonium import output
from quarkonium import plotter_and_normaliser
from quarkonium_plotter import sch_solver
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
    (2,0): [0.5, 0, 1.3],
    (3,0):(1.1,0,2.0)
}


#need to update these real values with proper data
real_energy_c = {
    (1,0): 0.4380, 
    (1,1): 0.9,
    (1,2): 1.2,
    #(2,0): 1.1
}

energy_ranges_b = {
    (1,0): [1, 0, 1.5],
    (1,1): [1.5,0,1.6], 
    (1,2): [1.6, 0, 1.9],
    (2,0): [1.4, 0, 1.7],
    (2,1): [1.7,0,2.0],
    (3,0): [1.8, 0, 2.2],
    (3,1): [2.1, 0, 2.3],
    (4,0): [2.1, 0, 2.3]

}

#taking the ground states
real_energy_b = {
    (1,0):1.0327 , 
    (1,1):1.4934 ,
    (1,2): 1.7977,
    (2,0): 1.633,
    (2,1): 1.866,
    (3,0): 1.989,
    (3,1): 2.147,
    (4,0): 2.213,
}

energy_b_plot = [
    1.0327 , 
    1.633,
    1.989,
    #2.213,
]

energy_b = [
    1.4934 ,
    1.7977,
    1.866,
    2.147,
]


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

    
#charmonium(3,0)
#bottonium(4,0)

def normalise(u,r):
    integral = sp.integrate.simpson(u**2,r)                           
    normalised_u = u/(np.sqrt(integral))
    return normalised_u




hyperfine = True
bottonium_true = True
#making the plots
def plots():
    output_r = []
    output_u = []
    #fig, axs = plt.subplots(ncols = 2, nrows = 2, figsize=(10, 4))
    

    fig = plt.figure(figsize = (10,4))

    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[6, 2],
        height_ratios=[1, 1]
    )

    ax1 = fig.add_subplot(gs[0, 0])   # top-left
    ax2 = fig.add_subplot(gs[1, 0])   # bottom-left
    ax3 = fig.add_subplot(gs[:, 1])   # right column spanning both rows

    plt.subplots_adjust(wspace=0, hspace=0.03)
    
    
    names = ['1S', '1P', '1D']
    colours = ['firebrick', 'chocolate', 'gold']
    colours = ['#81171B', '#C75146', '#EA8C55']
    colours = ['#E9C46A', '#C75146', '#81171B']
    colours = ["#EDC45F", 'chocolate', "#911A1E"]
    colours = ["#FFB700", "#D9673A", '#81171B']
   
    names_b = ['1S', '2S', '3S']
    #colours_b = ['royalblue', 'cornflowerblue', 'lightblue']
    colours_b = ['firebrick', 'cornflowerblue', 'mediumseagreen']
    colours_b = ['firebrick', "#7471B1", "#599AB8" ]
    colours_b = ['#EDC45F', "#36A8BA", '#264653']
    colours_b = ["#FFB700", "#36A8BA", '#264653']
    #colours_b = ['#81171B', '#006868',]
    sizes = [2.5,1,1]

    #analytic_c = [0.4381, 0.86871, 1.227]
    #analytic_b = [1.0327,  1.633, 1.989]
    #for E_c, E_b in zip(analytic_c, analytic_b):
    #    ax3.hlines(y=E_c, xmin=0.075, xmax=0.425, linewidth = 1, color = 'black')
    #    ax3.hlines(y=E_b, xmin=0.575, xmax=0.92, linewidth = 1, color = 'black')

    for n_l, name, colour, size in zip(energy_ranges_c, names, colours,sizes):
        E, u, r = charmonium(n_l[0], n_l[1])

        ax3.hlines(y=E, xmin=0.075, xmax=0.425, linewidth = 2, color = colour)

    
        #axs[0,1].axvspan(0.15, 0.85, ymin = E - 0.05*E, ymax = E + 0.05*E, color=colour, alpha=0.5)
        ax3.axhspan(E - 0.05*E, E + 0.05*E, xmin=0.075, xmax=0.425, color='lightgrey', alpha = 0.5,edgecolor = 'none')
        #print(f'this colour {colour} should have range',  E - 0.05*E,  E + 0.05*E)
        
        ax3.text(0.18, E+0.06, f'{name}', va='center', weight = 450, fontsize = 12)
        ax3.set_ylabel('Binding energy $E_{nl}$', labelpad = 12)
        #axs[0,1].set_ylim(0.3, 1.4)
        #ax3.set_ylim()

        ax1.scatter(r, u, label = str(n_l), marker = '.', s = size, color = colour)
        
        print(n_l,'was plotted successfully')
        output_r.append(r)
        output_u.append(u)
        #axs[0].set_aspect(6)  
        #axs[0,0].legend()
        #axs[0,0].set_xlabel('Separation $r$ (GeV)')
        ax1.set_ylabel('$u_{nl}(r)$', labelpad = 10)
        #axs[0,0].set_xticks([])
    j_psi = 0.437 + 0.160
    ax3.hlines(y=0.437 + 0.160, xmin=0.075, xmax=0.425, linewidth = 2, color = colours[0], linestyle =(0, (5, 1)))
    ax3.axhspan(j_psi - 0.05*j_psi, j_psi + 0.05*j_psi, xmin=0.075, xmax=0.425, color='lightgrey', alpha = 0.5, edgecolor = 'none')
    ax3.text(0.18, j_psi+0.06, '1S', va='center', weight = 450, fontsize = 12)
    #axs[0,1].hlines(y=1.20, xmin=0.3, xmax=0.7, colors='blue')
    ax3.set_xticks([])
    ax3.yaxis.tick_right()
    ax3.yaxis.set_label_position("right")
    ax3.set_xlim(0, 1)
    

    hyperfine_splittings = [0.047178795, 0.03892257,0.0367589]
    if bottonium_true:
        for E, name, colour,size, hyper in zip(energy_b_plot, names_b, colours_b, sizes,hyperfine_splittings):
            ax3.hlines(y=E, xmin=0.575, xmax=0.92, linewidth = 2, color = colour)
            ax3.text(0.5 + 0.18, E+0.10, f'{name}', va='center', weight = 450, fontsize = 12)
            ax3.set_ylabel('Binding energy $E_{nl}$')
            ax3.axhspan(E - 0.05*E, E + 0.05*E, xmin=0.575, xmax=0.92, color='lightgrey', alpha=0.5, edgecolor = 'none')
            
            ax3.hlines(y=E+hyper, xmin=0.575, xmax=0.92, linewidth = 2, linestyle = '--', color = colour)
            #ax3.axhspan(E+hyper - 0.05*E+hyper, E+hyper + 0.05*E+hyper, xmin=0.575, xmax=0.92, color='lightgrey', alpha=0.5, edgecolor = 'none')
       # ax3.set_ylim(0.95, 2.2)

        _, _, u_1, r_1 = sch_solver(0,4.183,4.183,0.33, 1.348327636718755, 1.0327, 2.7) #n = 1
        _, _, u_2, r_2 = sch_solver(0,4.183, 4.183, 0.33, 0.71376953125, 1.633, 5.3) #n = 2 
        _, _, u_3, r_3 = sch_solver(0,4.183, 4.183, 0.33, 0.71376953125, 2.42474365234375, 5.5) #n = 3 
        #_, _, u_4, r_4 = sch_solver(0,4.183, 4.183, 0.33, 0.71376953125, 3.09832763671875, 6.3) #n = 4

        ax2.hlines(y=0, xmin = min(r), xmax = 11.5, linestyle = '--', color = 'grey')
        ax2.scatter(r_1, u_1, marker = '.', s = sizes[0], label = 'n = 1', color = colours_b[0])
        ax2.scatter(r_2, u_2, marker = '.', s = sizes[1], label = 'n = 2', color = colours_b[1])
        ax2.scatter(r_3, u_3, marker = '.', s = sizes[2], label = 'n = 3', color = colours_b[2])
        #axs[1].scatter(r_4, u_4, marker = '.', s = 2, label = 'n = 4')
        #axs[1].set_aspect(2)  
        #axs[1,0].legend()



        ax2.set_ylabel('$u_{nl}(r)$', labelpad = 3)
        ax2.set_xlabel('Separation $r$ (GeV)')
        

    #axs[1,1].yaxis.tick_right()
    #axs[1,1].yaxis.set_label_position("right")
    #axs[1,1].set_xticks([])
    #axs[1,1].set_xlim(0, 1)
    ax1.tick_params( direction = 'in', width = 1.4, length = 6)
    ax2.tick_params( direction = 'in', width = 1.4, length = 6)
    ax1.minorticks_on()
    ax2.minorticks_on()
    ax1.tick_params( which =  'minor', direction = 'in') 
    ax2.tick_params( which =  'minor', direction = 'in') 
    ax1.set_xlim(xmin = 0, xmax = 11.5)
    ax1.set_ylim(ymin = 0)
    ax2.set_xlim(xmin = 0, xmax = 11.5)


    ax1.text(11.3, 0.65, r"$c\bar{c}$",ha='right', va='top', fontsize = 12, weight = 400)
    ax2.text(11.3, 1.05, r"$b\bar{b}$",ha='right', va='top', fontsize = 12, weight = 400)

    plt.minorticks_on()
    plt.tick_params( direction = 'in', width = 1.4, length = 6)
    plt.tick_params( which =  'minor', direction = 'in') 

    plt.savefig('summative/results_test_minor.png', bbox_inches = 'tight')
    plt.savefig('summative/results_test_minor.svg', bbox_inches = 'tight')
    #plt.show()
    


#Hyperfine splitting for n = 1 is 0.1604870964121239     
#Hyperfine splitting for n = 2 is 0.11566129414755096
plots()


 #not interesting to do hyperfine, radial wavefunction does not look any good
'''if hyperfine: 
        hf_1 = 0.1604870964121239  
        _ , _, _, _, _, final_node = sch_solver(0, c.m_c, c.m_c, E + hf_1 , c.alpha_c, c.beta_c, rmax)
        nodes_nb, turning_points_nb, u_1, v, r_1, final_node = sch_solver(0, c.m_c, c.m_c, E + hf_1, c.alpha_c, c.beta_c, final_node)
        integral = sp.integrate.simpson(u_1**2,r_1)                           
        normalised_u = u_1/(np.sqrt(integral))
        #normalised_check = sp.integrate.simpson(normalised_u**2, r_1)
        axs[0].scatter(r_1, normalised_u)'''
'''
    for n_l in energy_ranges_b:
        E, u, r = bottonium(n_l[0], n_l[1])
        axs[0].scatter(r, u, label = str(n_l), marker = '.', s = 2)
        print(n_l,'was plotted successfully')
        output_r.append(r)
        output_u.append(u) 
        axs[1].set_aspect(6)
        axs[1].legend()
        axs[1].set_xlabel('Separation $r$ (GeV)')
        axs[1].set_ylabel('$u_{nl}(r)$')

    '''


