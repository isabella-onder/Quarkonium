from quarkonium import output
from quarkonium import sch_solver
import cextension as ce
import bextension as be
import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import cmath
import machine 
import constants as c
import random

#here I want to make the code to be able to make a spectrum: that is, given a starting population and the various modes of decay
#for now: copy pasting the homework from last time

#motivate it as seeing how much of what and how much j/psi in the context of suppression QGP (can do one with real values and one with mine)

def has_transitioned(prob):
    r = random.random()
    if r < prob:
        return True
    else:
        return False


def evolveOne(currentState, rules):
    # collect all rules that apply
    applicable = [(final, prob) for (initial, final, prob) in rules
                  if initial == currentState]

    if not applicable:
        return currentState

    r = random.random()
    cumulative = 0.0

    for final, prob in applicable:
        cumulative += prob
        if r < cumulative:
            return final

    return currentState


def evolveMany(states, rules):
    newState = []
    for k in range (len(states)):
        newState.append(evolveOne(states[k], rules))
    return newState

def evolve_system(NA, NB, NC, ND, NE, rules, n_steps):
    
    state = (['A'] * NA)+(['B'] * NB)+(['C'] * NC)+(['D'] * ND)+(['E'] * NE)
    
    A_count = []
    B_count = []
    C_count = []
    D_count = []
    E_count = []
    
    for i in range (n_steps +1):
        a = 0
        b = 0
        c = 0
        d = 0
        e = 0
        for k in range(len(state)):
            if state[k] == 'A':
                a = a + 1
            elif state[k] == 'B':
                b = b + 1
            elif state[k] == 'C':
                c = c +1
            elif state[k] == 'D':
                d = d +1
            elif state[k] == 'E':
                e = e +1
            else:
                print('that seems not to be a valid input')
        A_count.append(a)
        B_count.append(b)
        C_count.append(c)
        D_count.append(d)
        E_count.append(e)
        state = evolveMany(state,rules)

    # YOUR CODE HERE
    return np.array(A_count), np.array(B_count), np.array(C_count), np.array(D_count), np.array(E_count)

##############################################################################################################
n_steps = 200
t_total = 0.00000000001
dt = t_total/n_steps

xs = []
for k in range((n_steps)+1):
    xs.append(dt*k)

xs1 = []
for k in range((n_steps)):
    xs.append(t_total + dt*(k+1))



##############################################################################################################
#here I will put in all the half lifetimes##########################################################################
##############################################################################################################
sec = 6.5821195695091*10**(-16+9) #converting eV to seconds
width = 92.6 *10**3 #the total width in eV

#inputting percentages from pdg live - by def of a width, is already related to lambda
t_av_hadrons = 1/(0.1346 * width) * sec 
print(t_av_hadrons)
t_av_threeg = 1/(0.641 * width) * sec
t_av_onetwo = 1/(0.088 * width) * sec
t_av_positron = 1/(0.05971 * width) * sec
t_av_muon = 1/(0.05961 * width) * sec

#t_av_mag = 1/(0.0141 * width) * sec


p_B = 1 - np.exp(-dt/t_av_hadrons)
p_C = 1 - np.exp(-dt/t_av_threeg)
p_D = 1 - np.exp(-dt/t_av_onetwo)
p_E = 1 - np.exp(-dt/t_av_positron)
p_F = 1 - np.exp(-dt/t_av_muon)


##############################################################################################################
#here are all the rules: who decays into what with which probability (just need to match and assign them) ############################################
##############################################################################################################
rules_1 = [
    ('A','B',p_B),
    ('A','C',p_C),
    ('A','D',p_D),
    ('A','E',p_E),
    #('C','A', p_C)   
]

def run() : 
    ##############################################################################################################
    #here i will put in all the possible populations - they all start with 0 except for J/psi
    ##############################################################################################################
    NA = 1000
    NB = 0
    NC = 0
    ND = 0
    NE = 0
    y_A, y_B, y_C, y_D, y_E = evolve_system(NA, NB, NC, ND, NE, rules_1, 2* n_steps)

    
    return y_A,y_B,y_C,y_D,y_E

y_a_tot,y_b_tot,y_c_tot, y_d_tot, y_e_tot = run()
plt.figure(figsize=(8, 4))
#plt.ylim(0,270)
plt.xlabel('time (h)')
plt.ylabel('Number of particles in given state')
plt.plot(xs+xs1, y_a_tot,label = 'A' )
plt.plot(xs+xs1, y_b_tot,label = 'B' )
plt.plot(xs+xs1, y_c_tot,label = 'C' )
plt.plot(xs+xs1, y_d_tot,label = 'D' )
plt.plot(xs+xs1, y_e_tot,label = 'E' )
plt.legend(loc = 'center right')
plt.title('Evolution of particle count according to state')
plt.show()