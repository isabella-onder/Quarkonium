from quarkonium import energy_range_finder
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import machine 
import constants as c

#for the shiggles: spherical harmonics calculator, to be able to produce full wavefunction
def spherical_harmonics(theta, phi, l, m):
    legendre_pol = sp.special.legendre(np.cos(theta))
    y_lm = (-1)**m * ((2*l+1)/(4 * np.pi) * np.math.factorial(l-m)/np.math.factorial(l+m))**(1/2) * legendre_pol * cmath.exp(complex(0,m*phi))
    return y_lm

#though we will only be studying l = 0 states
y_00 = np.sqrt(1/(4*np.pi))

#finding the wavefunction at the origin for CHARMONIUM. Assume that necessarily l = 0
def origin_c(n):
    _, _, u, r = energy_range_finder(0, c.m_c, c.m_c, machine.get_energy_range(n, 0), c.alpha_c, c.beta, c.rmax)
    print('this is the first 2 u', u[0], u[1])
    print('with corresponding 2 r', r[0], r[1])



origin_c(1)