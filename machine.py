from quarkonium import energy_range_finder
from quarkonium import plotter_and_normaliser

#assume alpha is fixed and beta is already found from the beta file
alpha = 0.40
beta = 0.210083

m_c = 1.2730 #Gev.c^-2

rmax = 30 #because it has been sufficient for now

#explicitely write and try energy range dictionary - arbitrary ranges
#(could make it like my hydrogen though that would assume known peaks)
energy_ranges = {
    (1,0): [0,0,0.5], 
    (1,1): [0,0,1],
    (1,2): [0.8, 0, 1.5],
    (2,0): [0.5, 0, 1]
}

#need to update these real values with proper data
real_energy= {
    (1,0): 0.4380, 
    (1,1): 0.9,
    (1,2): 1.2,
    (2,0): 1.1,
}


def get_energy_range(n, l):
    # Lookup the energy range for the given n and l
    try:
        return energy_ranges[(n, l)]
    except KeyError:
        raise ValueError(f"No energy range found for n={n}, l={l}. Please check the input.")
    

def charmonium(n,l):
    return energy_range_finder(l,m_c,m_c,get_energy_range(n,l), alpha, beta, rmax)


    
charmonium(1,0)