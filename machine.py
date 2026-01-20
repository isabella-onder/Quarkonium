from quarkonium import energy_range_finder
from quarkonium import output
from quarkonium import plotter_and_normaliser
import constants as c


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
    (2,0): [1.4, 0, 1.7]


}

#need to update these real values with proper data
real_energy_c = {
    (1,0): 0.4380, 
    (1,1): 0.9,
    (1,2): 1.2,
    (2,0): 1.1,
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
    E_n, _ , u, v, r =  output(l,c.m_c,c.m_c,get_energy_range(n,l,quark), c.alpha_c, c.beta, rmax)  #returns the final E_n and corresponding final_node (i.e. where it was cut off I think)

    print(f'\n##############################################################################\n {quark} QUARKONIUM \n The final estimated energy for n = {n}, l = {l} is E_nl = ', E_n, 'GeV', '\n############################################################################## \n' )

def bottonium(n,l):
    quark = 'BOTTOM'

    E_n, _ , u , v, r = output(l, c.m_b, c.m_b, get_energy_range(n,l,quark), c.alpha_b, c.beta, rmax )
    print(f'\n##############################################################################\n {quark} QUARKONIUM \n The final estimated energy for n = {n}, l = {l}  is E_nl = ', E_n, 'GeV', '\n############################################################################## \n' )


    
charmonium(1,0)
#bottonium(2,0)