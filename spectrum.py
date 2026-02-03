from quarkonium import output
from quarkonium import sch_solver
import extension as e
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

'''
def evolveOne(currentState, rules):

#this is to ensure that if it is not in rules, will remain unchanged
    comparatives = [this[0] for this in rules]
    i = 0
    for k in range (len(comparatives)):
        if currentState == comparatives[k]:
            i = i + 1
    if i == 0:
        return currentState
    
#this is to make it change according to which rule applies to it
    for i in range(len(rules)):
        initial, final, prob = rules[i]
        if currentState == initial:
            outcome = has_transitioned(prob)
            if outcome == True:
                currentState = final
                return currentState
            else:
                return currentState
        else:
            currentState = currentState
'''

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

def evolve_system(NA, NB, NC, rules, n_steps):
    
    state = (['A'] * NA)+(['B'] * NB)+(['C'] * NC)
    
    A_count = []
    B_count = []
    C_count = []
    
    for i in range (n_steps +1):
        a = 0
        b = 0
        c = 0

        for k in range(len(state)):
            if state[k] == 'A':
                a = a + 1
            elif state[k] == 'B':
                b = b + 1
            elif state[k] == 'C':
                c = c +1
                #print('we have indeed added a state to c')
            else:
                print('that seems not to be a valid input')
        A_count.append(a)
        B_count.append(b)
        C_count.append(c)
        state = evolveMany(state,rules)

    # YOUR CODE HERE
    return np.array(A_count), np.array(B_count), np.array(C_count)

##############################################################################################################
n_steps = 200
t_total = 100
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

#inputting percentages from pdg
t_half_hadrons = 0.1346
t_half_threeg = 15.7
t_half_C = 3.2

#need to convert to averages
t_av_A = 100.1/np.log(2)
t_av_B = 100.7/np.log(2)
t_av_C = 100.2/np.log(2)

p_A = 1 - np.exp(-dt/t_av_A)
p_B = 1 - np.exp(-dt/t_av_B)
p_C = 1 - np.exp(-dt/t_av_C)

##############################################################################################################
#here are all the rules: who decays into what with which probability (just need to match and assign them) ############################################
##############################################################################################################
rules_1 = [
    ('A','C',p_B), 
    ('A','B',p_A),
]

def run() : 
    ##############################################################################################################
    #here i will put in all the possible populations - they all start with 0 except for J/psi
    ##############################################################################################################
    NA = 250
    NB = 0
    NC = 0
    y_A, y_B, y_C = evolve_system(NA, NB, NC, rules_1, 2* n_steps)

    
    return y_A,y_B,y_C

y_a_tot,y_b_tot,y_c_tot = run()
plt.figure(figsize=(8, 4))
#plt.ylim(0,270)
plt.xlabel('time (h)')
plt.ylabel('Number of particles in given state')
plt.plot(xs+xs1, y_a_tot,label = 'A' )
plt.plot(xs+xs1, y_b_tot,label = 'B' )
plt.plot(xs+xs1, y_c_tot,label = 'C' )
plt.legend(loc = 'center right')
plt.title('Evolution of particle count according to state')
plt.show()