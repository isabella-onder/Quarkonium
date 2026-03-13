from quarkonium import output
from quarkonium import sch_solver
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import machine 
import constants as c

plot_wanted = False 
quark = 'CHARM'

wrong_origin = False

#all the data below will initially be extracted in origin_c out of output from quarkonium.py and then reused independently from there

#for the shiggles: spherical harmonics calculator, to be able to produce full wavefunction
def spherical_harmonics(theta, phi, l, m):
    legendre_pol = sp.special.legendre(np.cos(theta))
    y_lm = (-1)**m * ((2*l+1)/(4 * np.pi) * np.math.factorial(l-m)/np.math.factorial(l+m))**(1/2) * legendre_pol * cmath.exp(complex(0,m*phi))
    return y_lm

#though we will only be studying l = 0 states
y_00 = np.sqrt(1/(4*np.pi))

#the widths as given by pdg: quoting without the errors, in GeV
#for j/psi specifically

total_width = 92.6 * 10**(-6)
three_gluons_percent = 0.641
one_g_two_g_percent = 0.088
lepton_positron_percent = 0.05971
muon_positron_percent = 0.05961
#those with 4 leptons are of the order of e-5 percent
hyperfine_percent = 0.0141
three_photons_percent = 1.16*10**(-5)  #very unlikely already as is, many other improbable things are more likely



#finding the wavefunction at the origin for CHARMONIUM. 
# Assume that necessarily l = 0 
def origin_c(n):
    type_c = 'CHARM'
    _, _, u, v, r = output(0, c.m_c, c.m_c, machine.get_energy_range(n, 0, type_c), c.alpha_c, c.beta_c, c.rmax)
    origin = u[0], v[0]
    #just_after_origin = u[1], v[1]
    #print('with corresponding 2 r', r[0], r[1], r[5], r[1000])
    #print('this is the first 2 u', u[0], u[1], u[5], u[10000])
    print('this is the first 2 v ', v[0], v[1], v[5], v[1000])
    return origin, u, v, r

#origin_c(1)




#looking at the hyperfine splitting between spin up and down states, l = 0 necessarily
#use formula from booklet, approximating R(0) = u(0)/r = v(0)
def hyperfine_splitting(n):
    origin_values, _, _, _ = origin_c(n)
    u_0, v_0 = origin_values
    #v_0 = 0.7699520112846642
    delta_e = 8/9 * c.alpha_c/c.m_c**2 * v_0**2 
    print('these are the origin values u, v', origin_values)
    print('Hyperfine splitting for n = '+ str(n)+' is', delta_e)
    return delta_e

#hyperfine_splitting(2)

#inputting the formulae to calculate simple transitions: nothing more to compute
def extracting_mass(n = int, l = int, hyperfine = bool): #give n and l value, hyperfine says whether to add the hyperfine energy split
    E_n, _ , _, _, _ =  output(l,c.m_c,c.m_c,machine.get_energy_range(1, 0, quark), c.alpha_c, c.beta_c, c.rmax) #no care about rmax since just extracting the energy anyways
    if hyperfine:
        extracted_mass = E_n + hyperfine_splitting(1) + 2 * c.m_c
    else:
        extracted_mass = E_n + 2*c.m_c
    print('this is the calculated from scratch j_psi mass', extracted_mass, 'expect 3.1667GeV ')
    return extracted_mass

#at the start of each one I would do the same values called hence a function:
def header(n):
    M = extracting_mass(n,0,True)
    origin_values, _, _, _ = origin_c(n)
    energy, _, _ = machine.charmonium(n,0)
    
    u_0, v_0 = origin_values #we approx v_0 as R(0)
    #getting the correct v0
    _, _, _, _, _, lep_final_node_1 = sch_solver(0, c.m_c, c.m_c, energy + hyperfine_splitting(n) , c.alpha_c, c.beta_c, c.rmax)
    _, _, lep_u, lep_v, lep_r, lep_final_node = sch_solver(0, c.m_c, c.m_c, energy + hyperfine_splitting(n) , c.alpha_c, c.beta_c, lep_final_node_1)
    lep_integral_to_normalise = sp.integrate.simpson(lep_u**2,lep_r)                   
    lep_normalised_u = lep_u/(np.sqrt(lep_integral_to_normalise))
    print('these are the first values lep_v and hyperv', lep_v[0]/lep_integral_to_normalise)
    print('whereas this is the original v_0', v_0)
    lep_normalised_v = lep_v/lep_integral_to_normalise

    print(M, energy, v_0, lep_normalised_v[0])
    return M, energy, v_0, lep_normalised_v[0]


#header(2)

#because we are only doing the same, J/psi every time: no need to compute every time. 
# If necessary, change header and then change these values
M = 3.1434377800058737 #GeV total mass of the state
energy = 0.43695068359375 #Gev binding energy of ground state as found by solver
v_0 = 0.9148623504190927 #value at origin for eta_c(1S)
lep_normalised_v = 1.3096702583455615 #value at origin for J/psi
#v_0 = 0.7442660513414027 #value at origin for eta_c(2S)
#lep_normalised_v = 
#testing to see original values: taking M = 2*m_c and v(0) for eta_c as wavefunction at origin
radial_at_origin_eta = 0.7699520112846642 #value at origin for eta_c

#for S-2 instead
M = 3.143051527611561
energy = 1.09658203125
v_0 = 0.7442660513414027
lep_normalised_v = 0.8498642794425497

if wrong_origin:
    M = 2 * c.m_c
    lep_normalised_v = radial_at_origin_eta

#THE FOLLOWING ARE ALL assuming for J/psi specifically, so n = 1, l = 0
def lepton_decay(alpha):
    Psi = lep_normalised_v * y_00
    #M, energy, v_0, lep_normalised_v = header()
    lepton_width_o = 4 * (1/137)**2 * c.e_c**2/M**2 * v_0**2 *(1 - 16/3 * alpha/np.pi) #using the full version from one of the papers (though still not entire)
    lepton_width = 4 * (1/137)**2 * c.e_c**2/M**2 * lep_normalised_v**2 *(1 - 16/3 * alpha/np.pi) #using the full version from one of the papers (though still not entire)
    alternative_width = (16*np.pi *(1/137)**2*c.e_c**2)/M**2 * (1-4*c.m_c**2/M**2)**(1/2) * (1+2*c.m_c**2/M**2) * Psi**2
    print('this is the lepton width', lepton_width,  'whereas the real one is', lepton_positron_percent*total_width)
    print('this is the alternative lepton width', alternative_width)
    return lepton_width
#lepton_decay(c.alpha_c)


def three_gluons(alpha):
    
    
    #the simplification one
    three_gluon_width_simple_o = 40*(np.pi**2 - 9)/(81*np.pi) * alpha**3/M**2 * lep_normalised_v**2*(1-3.7*alpha/np.pi) 
    three_gluon_width_simple = 40*(np.pi**2 - 9)/(81*np.pi) * alpha**3/M**2 * lep_normalised_v**2*(1-3.7*alpha/np.pi) 
    #* (1 - 3.7 * c.alpha_c/np.pi)
    print('this is the three gluon width', three_gluon_width_simple)

    #the full one from the paper with its values, putting all in GeV
    #includes relativistic and radiative corrections, which above does not
    N = 0.57 * 10**(3/2) #from the paper
    beta = 0.310 #from the paper
    kappa = 3*(112+25*np.pi**2)/(16*(np.pi**2-9)) #as defined in the paper
    #print('should find that this is approx 0.77', kappa*beta**2/M**2)
    C = 3.7 #from the paper

    three_gluon_width = (80*(np.pi**2-9)*alpha**3 * N**2 * beta**3)/(81*np.pi**(9/2)*M)*(1-kappa*beta**2/M**2)*(1-C*alpha/np.pi)
    print('this is the paper gluon width', three_gluon_width)
    print('this is the experimental gluon width', three_gluons_percent * total_width)

    return three_gluon_width

#three_gluons(c.alpha_c)

def one_g_two_g(alpha): #i.e. one photon (gamma) and two gluons


    one_two_width_o = 32*(np.pi**2-9)/9 *c.e_c**2*(1/137)*c.alpha_c**2/c.m_c**2*v_0**2
    one_two_width = (32*(np.pi**2-9)*alpha**2*(1/137))/(9*np.pi*M**2) * lep_normalised_v**2 * (1-6.7*alpha/np.pi)
    print('this is the width for one photon and two gluons', one_two_width)
    print('this is the original width for one photon and two gluons', one_two_width_o)

    print('this is the experimental width for one photon two gluons', one_g_two_g_percent * total_width)
    return one_two_width

#one_g_two_g(c.alpha_c)

def three_photons():

    

    three_photon_width_o = 8*(np.pi**2-9)/(9) * c.e_c**2*(1/137)**3/M**2 * v_0**2 #I clearly just added a /pi 
    three_photon_width_c = 8*(np.pi**2-9)/(9) * c.e_c**2*(1/137)**3/M**2 * lep_normalised_v**2 #I clearly just added a /pi 
    print('this is the three photon width original', three_photon_width_o)
    print('this is the three photon width correct', three_photon_width_c)
    print('this is the real experimental value', total_width * three_photons_percent)
    return three_photon_width_c

#three_photons()



#this is not working very well: let's look at ratio formulae instead (need to fact check those formulae)


#transitioning magnetically width: use formula given in paper
#this works only for when the transitioning n are equal: otherwise, cannot extract photon momentum
def mag_transition(n):
    type_c = 'CHARM'
    #constants for charm quark
    alpha_em = 1/137
    e_Q = 2/3
    kappa_Q = 4/3 * c.alpha_c/(2*np.pi)


    k_gamma = (hyperfine_splitting(n)) #assuming it is indeed in N.U., otherwise need to divide by c
    k_gamma = 0.11
    
    #INTEGRAL from 0 to inf of R_n and R_n' * r^2:
    integral = 1 #since our n are equal, and u normalised, this is 1



    gamma_width = 4 * alpha_em * e_Q**2 /(3*(c.m_c)**2) * (1+kappa_Q)**2 *k_gamma**3
    print('Magnetic transition width', gamma_width, 'in GeV')
    return gamma_width

#mag_transition(1) #indeed, it overestimates it!!!!!!!! 


#doing that integral for the hyperfine splitting, as in the booklet: int(RR'r**2 dr). It seems like from notation that those R
#only depend on n', but since step size the same, can just multiply the u and sum over. 
def hyper_integral(u_lower, u_higher, r_lower, r_higher):
    mixed_u = [u_1 *u_2 for u_1, u_2 in zip(u_lower, u_higher)]
    #mixed_r = [r_1 * r_2 for r_1, r_2 in zip(r_higher, r_lower)]
    
    
    if np.array_equal(r_lower, r_higher):
        print('the r arrays are indeed the same')
    else:
        print('the r arrays are not the same, do not proceed to panic')
        print(r_lower[-1], r_higher[-1])
    

    integral_putain = sp.integrate.simpson(mixed_u, r_lower)
    print('this is da integral', integral_putain)
    return integral_putain


#the same as mag_transition but accounting for the fact that the energy has changed in one, and therefore perhaps also 
#the values of its radial wavefunction. So all is the same except that the integral may differ (and we have to extract that splitting
#solution)
def alternative_mag_transition(n):
    type_c = 'CHARM'
    #constants for charm quark
    alpha_em = 1/137
    e_Q = 2/3
    kappa_Q = 4/3 * c.alpha_c/(2*np.pi)

    #need to input the energy here
    energy, _, _ = machine.charmonium(1,0) #GeV, as pdt by charmonium(1,0) in machine.py

    k_gamma = (hyperfine_splitting(n)) #assuming it is indeed in N.U., otherwise need to divide by c
    k_gamma = 0.11

  
    #hyperfine results: use l = 0 and other constants as usual. do it twice to get the final_node and then normalise
    _, _, _, _, _, hyper_final_node_1 = sch_solver(0, c.m_c, c.m_c, energy + k_gamma, c.alpha_c, c.beta_c, c.rmax)
    hyper_nodes_nb, hyper_turning_points_nb, hyper_u, hyper_v, hyper_r, hyper_final_node = sch_solver(0, c.m_c, c.m_c, energy + k_gamma, c.alpha_c, c.beta_c, hyper_final_node_1)
    hyper_integral_to_normalise = sp.integrate.simpson(hyper_u**2,hyper_r)                   
    hyper_normalised_u = hyper_u/(np.sqrt(hyper_integral_to_normalise))
    
    #usual result, except that we make it terminate at hyper final node such that we are summing over the same things
    #hence have to do it manually, for there not to be an issue with finding the energy
    #_, _, u, v, r = output(0, c.m_c, c.m_c, machine.get_energy_range(n, 0, type_c), c.alpha_c, c.beta, hyper_final_node)

    _, _, casual_u, casual_v, casual_r, casual_final_node = sch_solver(0, c.m_c, c.m_c, energy , c.alpha_c, c.beta_c, hyper_final_node_1)
    casual_integral_to_normalise = sp.integrate.simpson(casual_u**2,casual_r)                   
    casual_normalised_u = casual_u/(np.sqrt(casual_integral_to_normalise))
    print('these are the first values v and hyperv', casual_v[0]/casual_integral_to_normalise, hyper_v[0]/hyper_integral_to_normalise)

    if plot_wanted:
        fig, axs = plt.subplots(ncols = 2, nrows = 1)
        axs[0].scatter(casual_r, casual_u, marker = '.')
        axs[1].scatter(hyper_r, hyper_u, marker = '.')
        plt.show()

    integral = hyper_integral(casual_normalised_u, hyper_normalised_u, casual_r, hyper_r)
    
    #using the given formula for magnetic transition, using the previously found splitting energy (fine or hyperfine?)
    gamma_width = 4 * alpha_em * e_Q**2 /(3*(c.m_c)**2) * (1+kappa_Q)**2 *k_gamma**3 * integral
    print('Magnetic transition width', gamma_width, 'in GeV using the alternative mag transition')
    print('Magnetic transition width', gamma_width/total_width, 'in fraction form using the alternative mag transition')

    return gamma_width
#alternative_mag_transition(2)


#using the ratio formulae in a paper to see whether they are better approx
def ratios(alpha):
    three_g_ratio = (5*alpha**3)/(54*(1/137)**3*c.e_c**6) * (three_photons_percent*total_width)
    one_gamma_two_g_ratio = (2*alpha**2)/(3*(1/137)**2*c.e_c**4) * (three_photons_percent*total_width)

    three_g_ratio_e = (5*alpha**3)/(54*(1/137)**3*c.e_c**6) * three_photons()
    one_gamma_two_g_ratio_e = (2*alpha**2)/(3*(1/137)**2*c.e_c**4) * three_photons()


    print('these are the exact experimental three gluons & one photon and two gluons widths', total_width * three_gluons_percent, total_width * one_g_two_g_percent)
    print('these are my values three gluons & one photon and two gluons widths', three_g_ratio_e, one_gamma_two_g_ratio_e)
    print('these are the values found via ratio', three_g_ratio, one_gamma_two_g_ratio)
    ratio_list = [three_g_ratio, one_gamma_two_g_ratio]
    return ratio_list
#ratios(c.alpha_c)







