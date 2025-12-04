#then integrates and plots
import numpy as np
import matplotlib.pyplot as plt

from milestone import plot

mu = 0.000511 #in GeV-1
alpha = 1/137
a0 = 1/(mu*alpha)  # Bohr radius in GeV^-1
#a0 = 5.11 * 10**(-11) #bohr radius in m

colours = ['grey', '#B1D0ED', '#156082']

def analytic(r_num_array):

    

    #extracting the arrays used in the milestone code
    r_10 = r_num_array[0]
    r_20 = r_num_array[1]
    r_21 = r_num_array[2]


    
    
    
    #r = np.linspace(0.0, a0 * 10 ,1000) #hopefully also in GeV-1

    #using the analytic formulae from Sakurai appendix B
    #set Z = 1 for Hydrogen
    
    
   

    R_10 = (1/a0)**(3/2)*2*np.exp(-r_10/a0)
    u_10 = [r * R for r, R in zip(r_10, R_10)]
    u_10_squared = [u**2 for u in u_10]
 

    R_20 = (1/(2*a0))**(3/2)*(2-r_20/a0)*np.exp(-r_20/(2*a0))
    u_20 = [r * R for r, R in zip(r_20, R_20)]
    u_20_squared = [u**2 for u in u_20]

    R_21 = (1/(2*a0))**(3/2)*(r_21/(np.sqrt(3)*a0))*np.exp(-r_21/(2*a0))
    u_21 = [r * R for r, R in zip(r_21, R_21)]
    u_21_squared = [u**2 for u in u_21]


    

    fig, ax = plt.subplots()

    plt.plot(r_10/a0, u_10_squared, color = 'grey')
    plt.plot(r_20/a0, u_20_squared, color = '#B1D0ED')
    plt.plot(r_21/a0, u_21_squared, color = '#156082')
    #plt.savefig("analytic_plots.svg", bbox_inches = 'tight')
    #plt.show()

    #theoretical_array = np.array(u_10_squared, u_20_squared, u_21_squared)
    theoretical_array = [u_10, u_20, u_21]

    return theoretical_array
  
 


def comparator_with_residuals():


    #getting the numerical values u_nl and corresponding r arrays from milestone
    numerical_array, r_num_array = plot([[1,0], [2,0], [2,1]])
    #using the same r_arrays, getting the theoretical values (as calculated in Sakurai)
    theoretical_array = analytic(r_num_array)
    #

    

    fig1, axs_new = plt.subplots(3,1, sharex = True)
    print(axs_new.shape)
    fig2, ax_res = plt.subplots()
    for u_num,  u_theo, color, i, r_num in zip(numerical_array, theoretical_array, colours, [0,1,2], r_num_array):
        residual = [(num - theo)/theo for num, theo in zip(u_num, u_theo)]
        #print('these are residual', residual)
    

        u_num_squared = [u**2 for u in u_num]
        u_theo_squared = [u**2 for u in u_theo]


        #plotting: numerical, theoretical, residuals of u_nl
        axs_new[0].scatter(r_num/a0, u_num, color = color, marker = '.')
        axs_new[1].scatter(r_num/a0, u_theo, color = color, marker = '.')
        axs_new[2].scatter(r_num/a0, residual, color = color,  marker = '.', s = 1)

        #plotting blow up of residuals (at the start before it starts proper diverging)
        ax_res.scatter(r_num[10:int(len(r_num)/4)]/a0, residual[10:int(len(r_num)/4)], color = color, marker = '.', s = 1)
        #ax_res.set_ylim(np.min(residual), np.max(residual))

    #plt.show()
    axs_new[0].set_ylabel('Numerical')
    axs_new[1].set_ylabel('Analytic')
    axs_new[2].set_ylabel('Residuals \n ($u_{num} - u_{theo}$)')
    ax_res.set_ylabel('Residuals zoom \n ($u_{num} - u_{theo}$)')
    axs_new[2].set_xlabel('Separation r in $a_0$')
    fig1.savefig('comparator.svg', bbox_inches = 'tight')
    fig2.savefig('residual_zoom.svg', bbox_inches = 'tight')

comparator_with_residuals()
        