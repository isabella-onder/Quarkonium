from quarkonium import output
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import machine 
import constants as c

#all the data below will initially be extracted in origin_c out of output from quarkonium.py and then reused independently from there

#for the shiggles: spherical harmonics calculator, to be able to produce full wavefunction
def spherical_harmonics(theta, phi, l, m):
    legendre_pol = sp.special.legendre(np.cos(theta))
    y_lm = (-1)**m * ((2*l+1)/(4 * np.pi) * np.math.factorial(l-m)/np.math.factorial(l+m))**(1/2) * legendre_pol * cmath.exp(complex(0,m*phi))
    return y_lm

#though we will only be studying l = 0 states
y_00 = np.sqrt(1/(4*np.pi))

#finding the wavefunction at the origin for CHARMONIUM. 
# Assume that necessarily l = 0 
def origin_c(n):
    type_c = 'CHARM'
    _, _, u, v, r = output(0, c.m_c, c.m_c, machine.get_energy_range(n, 0, type_c), c.alpha_c, c.beta, c.rmax)
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
    delta_e = 8/9 * c.alpha_c/c.m_c**2 * v_0**2 
    print('these are the origin values u, v', origin_values)
    print('Hyperfine splitting for n = '+ str(n)+' is', delta_e)
    return delta_e

#hyperfine_splitting(1)

#transitioning magnetically width: use formula given in paper
#this works only for when the transitioning n are equal: otherwise, cannot extract photon momentum
def mag_transition(n):
    
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

mag_transition(1) #indeed, it overestimates it!!!!!!!! 