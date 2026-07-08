
import random
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from matplotlib.ticker import FuncFormatter



#here I want the code to just make the animation: keep it as a separate function of its own, as I may want to tweak things specifically

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

def evolve_system(NA, NB, NC, ND, rules, n_steps):
    
    state = (['A'] * NA)+(['B'] * NB)+(['C'] * NC)+(['D'] * ND)
    
    A_count = []
    B_count = []
    C_count = []
    D_count = []

    
    for i in range (n_steps +1):
        a = 0
        b = 0
        c = 0
        d = 0

        for k in range(len(state)):
            if state[k] == 'A':
                a = a + 1
            elif state[k] == 'B':
                b = b + 1
            elif state[k] == 'C':
                c = c +1
            elif state[k] == 'D':
                d = d +1
           
            else:
                print('that seems not to be a valid input')
        A_count.append(a)
        B_count.append(b)
        C_count.append(c)
        D_count.append(d)
        
        state = evolveMany(state,rules)

    # YOUR CODE HERE
    return np.array(A_count), np.array(B_count), np.array(C_count), np.array(D_count)

##############################################################################################################
n_steps = 200
t_total = 5
dt = t_total/n_steps

xs = []
for k in range((n_steps)+1):
    xs.append(dt*k)

xs1 = []
for k in range((n_steps)):
    xs.append(t_total + dt*(k+1))

#so that the plot does not have to do rounding
norm = 1e20

##############################################################################################################
#here I will put in all the half lifetimes##########################################################################
##############################################################################################################
sec = 6.5821195695091*10**(-25) #converting GeV to seconds


#inputting percentages from pdg live - by def of a width, is already related to lambda
t_av_onetwo = 1/(1.903609754614899e-06) * sec * norm
t_av_lepton = 1/(1.3650953620536091e-06) * sec * norm
t_av_threeg = 1/(4.5106026937836827e-05) * sec * norm


#B: one photon, two gluons
p_B = 1 - np.exp(-dt/t_av_onetwo)

#C: lepton
p_C = 1 - np.exp(-dt/t_av_lepton)

#D: three gluons
p_D = 1 - np.exp(-dt/t_av_threeg)



##############################################################################################################
#here are all the rules: who decays into what with which probability (just need to match and assign them) ############################################
##############################################################################################################
rules_1 = [
    ('A','B',p_B),
    ('A','C',p_C),
    ('A','D',p_D)

]

def run() : 
    ##############################################################################################################
    #here i will put in all the possible populations - they all start with 0 except for J/psi
    ##############################################################################################################
    NA = 1000
    NB = 0
    NC = 0
    ND = 0

    y_A, y_B, y_C, y_D = evolve_system(NA, NB, NC, ND, rules_1, 2* n_steps)

    
    return y_A,y_B,y_C,y_D

y_a_tot,y_b_tot,y_c_tot, y_d_tot = run()


fig, ax = plt.subplots()

plt.xlabel('Time ($10^{-19} s$)')
plt.ylabel('Percentage of particles in given state')


x_axis = xs+xs1
plot_A = ax.plot(x_axis, y_a_tot,label = 'ϒ(1S)', linewidth = 2.5, color = 'royalblue' )[0]
plot_B = ax.plot(xs+xs1, y_b_tot,label = 'γgg', linewidth = 2.5, color = 'sienna' )[0]
plot_C = ax.plot(xs+xs1, y_c_tot,label = 'lepton', linewidth = 2.5, color = 'goldenrod' )[0]
plot_D = ax.plot(xs+xs1, y_d_tot,label = 'ggg', linewidth = 2.5, color = 'firebrick' )[0]

#ax.legend(loc = ' right')
ax = plt.gca()
ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f'{x/10:g}'))
ax.legend()
plt.title('Evolution of particle count according to state')


def update(frame):
    # for each frame, update the data stored on each artist.
    x = x_axis[:frame]
    y = y_a_tot[:frame]

    # update the line plot:
    plot_A.set_xdata(x_axis[:frame])
    plot_A.set_ydata(y_a_tot[:frame])

    plot_B.set_xdata(x_axis[:frame])
    plot_B.set_ydata(y_b_tot[:frame])

    plot_C.set_xdata(x_axis[:frame])
    plot_C.set_ydata(y_c_tot[:frame])

    plot_D.set_xdata(x_axis[:frame])
    plot_D.set_ydata(y_d_tot[:frame])

    return (plot_A, plot_B, plot_C, plot_D)


ani = animation.FuncAnimation(fig=fig, func=update, frames=400, interval=30, repeat = False)

#if I want to save, just need to unhash
animation.Animation.save(ani, 'animations/formative_bottom_num.gif')
plt.show()

