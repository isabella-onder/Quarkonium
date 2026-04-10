from quarkonium import output
from quarkonium import sch_solver
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import machine 
import constants as c

plot_wanted = False 
quark = 'BOTTOM'

wrong_origin = False

#all the data below will initially be extracted in origin_b out of output from quarkonium.py and then reused independently from there

#for the shiggles: spherical harmonics calculator, to be able to produce full wavefunction
def spherical_harmonics(theta, phi, l, m):
    legendre_pol = sp.special.legendre(np.cos(theta))
    y_lm = (-1)**m * ((2*l+1)/(4 * np.pi) * np.math.factorial(l-m)/np.math.factorial(l+m))**(1/2) * legendre_pol * cmath.exp(complex(0,m*phi))
    return y_lm

#though we will only be studying l = 0 states
y_00 = np.sqrt(1/(4*np.pi))

#the widths as given by pdg: quoting without the errors, in GeV
#for j/psi specifically

total_width = 54.02 * 10**(-6)
tau_percent = 0.026
positron_percent = 0.0239
muon_percent = 0.0248
threegluons_percent = 0.817
onegtwog_percent = 0.022
jpsi_percent = 0.00054
eta_percent = 0.0294 #I am unsure this is the correct eta

threegluons_error = 0.007
positron_error = 0.0008
onegtwog_error = 0.006
eta_error = 0.0024



#finding the wavefunction at the origin for bottomONIUM. 
# Assume that necessarily l = 0 
def origin_b(n):
    type_b = 'BOTTOM'
    _, _, u, v, r = output(0, c.m_b, c.m_b, machine.get_energy_range(n, 0, type_b), c.alpha_b, c.beta_b, c.rmax)
    origin = u[0], v[0]
    #just_after_origin = u[1], v[1]
    #print('with corresponding 2 r', r[0], r[1], r[5], r[1000])
    #print('this is the first 2 u', u[0], u[1], u[5], u[10000])
    #print('this is the first 2 v ', v[0], v[1], v[5], v[1000])
    return origin, u, v, r

#origin_b(1)




#looking at the hyperfine splitting between spin up and down states, l = 0 necessarily
#use formula from booklet, approximating R(0) = u(0)/r = v(0)
def hyperfine_splitting(n):
    origin_values, _, _, _ = origin_b(n)
    u_0, v_0 = origin_values
    delta_e = 8/9 * c.alpha_b/c.m_b**2 * v_0**2 
    #print('these are the origin values u, v', origin_values)
    print('Hyperfine splitting for n = '+ str(n)+' is', delta_e)
    return delta_e

#hyperfine_splitting(4)

#inputting the formulae to calculate simple transitions: nothing more to compute
def extracting_mass(n = int, l = int, hyperfine = bool): #give n and l value, hyperfine says whether to add the hyperfine energy split
    E_n, _ , _, _, _ =  output(l,c.m_b,c.m_b,machine.get_energy_range(n, 0, quark), c.alpha_b, c.beta_b, c.rmax) #no care about rmax since just extracting the energy anyways
    if hyperfine:
        extracted_mass = E_n + hyperfine_splitting(1) + 2 * c.m_b
    else:
        extracted_mass = E_n + 2*c.m_b
    #print('this is the calculated from scratch bottom_1 mass', extracted_mass, 'expect 9.4604 GeV ')
    return extracted_mass
#extracting_mass(1,0,True)


#at the start of each one I would do the same values called hence a function:
def header(n):
    M = extracting_mass(n,0,True)
    origin_values, _, _, _ = origin_b(n)
    energy,_,_ = machine.bottonium(n,0)
    
    u_0, v_0 = origin_values #we approx v_0 as R(0)

    #getting the correct v0
    _, _, _, _, _, lep_final_node_1 = sch_solver(0, c.m_b, c.m_b, energy +hyperfine_splitting(n) , c.alpha_b, c.beta_b, c.rmax)
    _, _, lep_u, lep_v, lep_r, lep_final_node = sch_solver(0, c.m_b, c.m_b, energy +hyperfine_splitting(n), c.alpha_b, c.beta_b, lep_final_node_1)
    lep_integral_to_normalise = sp.integrate.simpson(lep_u**2,lep_r)                   
    lep_normalised_u = lep_u/(np.sqrt(lep_integral_to_normalise))
    print('these are the first values lep_v and hyperv', lep_v[0]/lep_integral_to_normalise)
    print('whereas this is the original v_0', v_0)
    lep_normalised_v = lep_v/lep_integral_to_normalise

    print(M, energy, v_0, lep_normalised_v[0])
    return M, energy, v_0, lep_normalised_v[0]

#header(2)

M = 9.43290539176862 #GeV total mass of B n = 1 spin up
energy = 1.03314208984375 #GeV energy of binding energy of n = 1 l = 0
v_0 = 1.540662693150416 # value of the wavefunction at the origin for n = 1 l = 0 spin down

#the one giving good results: keeping that one from the previous beta, since more appropriate alpha
lep_normalised_v = 2.8229989433130935 # value of the wavefunction at the origin n =  l = 0 spin up
#lep_normalised_v = 3.7371766755216824 

if wrong_origin:
    M = 2*c.m_b
    #lep_normalised_v = v_0

#THE FOLLOWING ARE ALL assuming for J/psi specifically, so n = 1, l = 0
def lepton_decay(alpha):
    
    Psi_o = v_0*y_00
    Psi = lep_normalised_v*y_00

    lepton_width = (16*np.pi*(1/137)**2*Psi**2/(9*(M)**2))*(1-(16*alpha)/(3*np.pi))
    lepton_width_original = (16*np.pi*(1/137)**2*Psi_o**2/(9*(M)**2))*(1-(16*alpha)/(3*np.pi))

    alternative_width = (16*np.pi *(1/137)**2*c.e_b**2)/M**2 * (1-4*c.m_b**2/M**2)**(1/2) * (1+2*c.m_b**2/M**2) * Psi**2
    
    print('############################################################## \n Lepton width original (GeV)', lepton_width_original, '\n #############################################################')
    print('############################################################## \n Lepton width correct (GeV)', lepton_width, '\n #############################################################')
    print('############################################################## \n Lepton width alternative (GeV)', alternative_width, '\n #############################################################')
    print('############################################################## \n Lepton width correct positron, tauons, muons (GeV)', positron_percent*total_width , tau_percent*total_width, muon_percent*total_width, '\n #############################################################')

    return alternative_width
#lepton_decay(0.17131313131313133)


def three_gluons(alpha):

    m = M
 
    Psi_o = v_0*y_00


    Psi = lep_normalised_v*y_00
    
    
    
    #the simplification one
    three_gluon_width_o = (160*(np.pi**2 - 9)*alpha**3*Psi_o**2)/(81*(m)**2)*(1-4.9*alpha/np.pi)
    three_gluon_width_correct = (160*(np.pi**2 - 9)*alpha**3*Psi**2)/(81*(m)**2)*(1-4.9*alpha/np.pi)
    print('############################################################## \n The three gluon width original (GeV)', three_gluon_width_o, '\n ##############################################################')
    print('############################################################## \n The three gluon width correct (GeV)', three_gluon_width_correct, '\n ##############################################################')
    print('############################################################## \n The three gluon width experimental (GeV)', threegluons_percent * total_width, '\n ##############################################################')

    return three_gluon_width_correct

#three_gluons()

def one_g_two_g(alpha): #i.e. one photon (gamma) and two gluons

    m = M
  
    Psi_o = v_0*y_00


    Psi = lep_normalised_v*y_00

    one_two_width_o = (128*(np.pi**2 - 9) * (1/137) * alpha**2 * Psi_o**2)/(81*(m)**2)*(1-1.7*alpha/np.pi)
    one_two_width_c = (128*(np.pi**2 - 9) * (1/137) * alpha**2 * Psi**2)/(81*(m)**2)*(1-1.7*alpha/np.pi)

    print('############################################################## \n One photon two gluon width (GeV)', one_two_width_o, '\n ##############################################################')
    print('############################################################## \n One photon two gluon width (GeV)', one_two_width_c, '\n ##############################################################')
    print('############################################################## \n One photon two gluon width (GeV)', onegtwog_percent * total_width, '\n ##############################################################')

    return one_two_width_c
#one_g_two_g()


#hmm cannot seem to see the experimental value: hold off for now
def three_photons():
    m = M
    Psi_o = v_0*y_00
    Psi = lep_normalised_v*y_00

    three_photons = (64*(np.pi**2 - 9) * alpha**3 * psi**2)/(2187*m**2)*(1 - 12.6 * alpha/np.pi)
    
    
    print('this is the three photon width', three_photons)
    print('this is the real experimental value', total_width * percent )
    return three_photons
#three_photons()





#doing the hyperfine splitting integral for the hyperfine splitting, as in the booklet: int(RR'r**2 dr). It seems like from notation that those R
#only depend on n', but since step size the same, can just multiply the u and sum over. 
def hyper_integral(u_lower, u_higher, r_lower, r_higher):
    mixed_u = [u_1 *u_2 for u_1, u_2 in zip(u_lower, u_higher)]
    #mixed_r = [r_1 * r_2 for r_1, r_2 in zip(r_higher, r_lower)]
    
    
    '''if np.array_equal(r_lower, r_higher):
        print('the r arrays are indeed the same')
    else:
        print('the r arrays are not the same, do not proceed to panic')
        print(r_lower[-1], r_higher[-1])
    '''

    integral_putain = sp.integrate.simpson(mixed_u, r_lower)
    #print('this is da integral', integral_putain)
    return integral_putain


#the same as mag_transition but accounting for the fact that the energy has changed in one, and therefore perhaps also 
#the values of its radial wavefunction. So all is the same except that the integral may differ (and we have to extract that splitting
#solution)
def alternative_mag_transition(n):
    type_b = 'bottom'
    #constants for bottom quark
    alpha_em = 1/137
    e_Q = -1/3
    kappa_Q = 4/3 * c.alpha_b/(2*np.pi)

    #need to input the energy here
    energy = 1.03314208984375 #GeV, as pdt by bottonium(1,0) in machine.py

    k_gamma = (hyperfine_splitting(n)) #assuming it is indeed in N.U., otherwise need to divide by c
    k_gamma = 0.11

  
    #hyperfine results: use l = 0 and other constants as usual. do it twice to get the final_node and then normalise
    _, _, _, _, _, hyper_final_node_1 = sch_solver(0, c.m_b, c.m_b, energy + k_gamma, c.alpha_b, c.beta_b, c.rmax)
    hyper_nodes_nb, hyper_turning_points_nb, hyper_u, hyper_v, hyper_r, hyper_final_node = sch_solver(0, c.m_b, c.m_b, energy + k_gamma, c.alpha_b, c.beta_b, hyper_final_node_1)
    hyper_integral_to_normalise = sp.integrate.simpson(hyper_u**2,hyper_r)                   
    hyper_normalised_u = hyper_u/(np.sqrt(hyper_integral_to_normalise))
    
    #usual result, except that we make it terminate at hyper final node such that we are summing over the same things
    #hence have to do it manually, for there not to be an issue with finding the energy
    #_, _, u, v, r = output(0, c.m_b, c.m_b, machine.get_energy_range(n, 0, type_b), c.alpha_b, c.beta, hyper_final_node)

    _, _, casual_u, casual_v, casual_r, casual_final_node = sch_solver(0, c.m_b, c.m_b, energy , c.alpha_b, c.beta_b, hyper_final_node_1)
    casual_integral_to_normalise = sp.integrate.simpson(casual_u**2,casual_r)                   
    casual_normalised_u = casual_u/(np.sqrt(casual_integral_to_normalise))
    #print('these are the first values v and hyperv', casual_v[0]/casual_integral_to_normalise, hyper_v[0]/hyper_integral_to_normalise)

    if plot_wanted:
        fig, axs = plt.subplots(ncols = 2, nrows = 1)
        axs[0].scatter(casual_r, casual_u, marker = '.')
        axs[1].scatter(hyper_r, hyper_u, marker = '.')
        plt.show()

    integral = hyper_integral(casual_normalised_u, hyper_normalised_u, casual_r, hyper_r)
    
    #using the given formula for magnetic transition, using the previously found splitting energy (fine or hyperfine?)
    gamma_width = 4 * alpha_em * e_Q**2 /(3*(c.m_b)**2) * (1+kappa_Q)**2 *k_gamma**3 * integral
    print('Magnetic transition width', gamma_width, 'in GeV using the alternative mag transition')
    return gamma_width
#alternative_mag_transition(1)


#using the ratio formulae in a paper to see whether they are better approx
def ratios(alpha):
    twogoneg_threeg_ratio = (4*(1/137))/(5*alpha)*(1-2.6*alpha/np.pi)
    threeg_muons_ratio = (10*(np.pi**2-9)*alpha**3)/(np.pi*(1/137)**2)*(M/(2*c.m_b)**2)**2*(1+0.43*alpha/np.pi)
    print('these are the ratio', twogoneg_threeg_ratio, threeg_muons_ratio)
    print('these are the experimental ratio', onegtwog_percent/threegluons_percent, threegluons_percent/muon_percent)
    ratio_list = [twogoneg_threeg_ratio, threeg_muons_ratio]
    return ratio_list
#ratios(c.alpha_b)







