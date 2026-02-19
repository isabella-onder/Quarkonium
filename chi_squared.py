import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import constants as c
import cextension as cex
import bextension as bex

from quarkonium import output
from quarkonium import sch_solver

#doing the chi squared

#choosing which mass to input
constituent_mass = False
if constituent_mass:
    m = 2*c.m_b #this is taking the mass as two quarks, i.e. constituenta
else:
    m = bex.extracting_mass(1,0,True) #this is getting the mass of the state including hyperfine

#taking those for which I have both a theoretical and experimental value, including ratios
#took the lepton decay as the positron percent arbitrarily, could have been muon



def psi():
    origin_values, _, _, _ = bex.origin_b(1)
    u_0, v_0 = origin_values #we approx v_0 as R(0)
    Psi = v_0*bex.y_00
    print(Psi)
    return Psi

Psi = 0.4346129216180489
Psi = 0.7963532990917851

#one photon two gluons
def one_g_two_g(alpha): #i.e. one photon (gamma) and two gluons


    one_two_width = (128*(np.pi**2 - 9) * (1/137) * alpha**2 * Psi**2)/(81*(m)**2)*(1-1.7*alpha/np.pi)
    return one_two_width


#positron
def lepton_decay(alpha):



    lepton_width = (16*np.pi*(1/137)**2*Psi**2/(9*(m)**2))*(1-(16*alpha)/(3*np.pi))
    return lepton_width

#three gluons
def three_gluons(alpha):

    

    three_gluon_width = (160*(np.pi**2 - 9)*alpha**3*Psi**2)/(81*(m)**2)*(1-4.9*alpha/np.pi)
    return three_gluon_width





#creating the alpha array, values along which to try
def values(type):
    start, end = 0, 0.8 #default values
    if type == 'CHARM':
        start = 0.1
        end = c.alpha_c + 0.1
    elif type == 'BOTTOM':
        start = 0
        end = c.alpha_b + 0.20

    step_nb = 100 # evert 0.01 is enough so do 30 steps for now to see trend
    alpha_arr = np.linspace(start, end,step_nb)
    print('This is the step size:', (end - start)/step_nb)
    return alpha_arr






def bex_chi_squared():
    bexp_widths = [ bex.onegtwog_percent*bex.total_width, bex.positron_percent*bex.total_width, bex.threegluons_percent*bex.total_width]
    btheo_widths_original = [bex.one_g_two_g(c.alpha_b), bex.lepton_decay(c.alpha_b), bex.three_gluons(c.alpha_b)]

    #calculate chi_squared for every given chi
    type_b = 'BOTTOM'
    alpha_arr = values(type_b)
    chi_squared_arr = []
    for alpha in alpha_arr:
        print('the next alpha is being computed')
        btheo_widths = [one_g_two_g(alpha),lepton_decay(alpha),three_gluons(alpha)]  #widths calculated for a given alpha
        if len(btheo_widths) != len(bexp_widths):
            print('There are not as many experimental as theo widths! need to update! ')
        #sum_elements = [(theo-exp)**2/error**2 for theo, exp, error in zip(btheo_widths, bexp_widths, bexp_widths_errors)] #
        sum_elements = [(theo-exp)**2 for theo, exp in zip(btheo_widths, bexp_widths)]
        chi_squared_value = sum(sum_elements)
        chi_squared_arr.append(chi_squared_value)

    idx_order = np.argsort(chi_squared_arr)
    arg_min_chi = idx_order[-1]
    #optimal_alpha = alpha_arr[arg_min_chi]
    optimal_alphas = [alpha_arr[i] for i in idx_order]
    optimal_alpha = optimal_alphas[0]

    btheo_widths_optimal = [one_g_two_g(optimal_alpha),lepton_decay(optimal_alpha),three_gluons(optimal_alpha)]
    print('###################################################################### \n','These are the optimal theo widths', btheo_widths_optimal, '\n #####################################################################')
    print('###################################################################### \n','These are the original theo widths', btheo_widths_original, '\n #####################################################################')
    print('###################################################################### \n','These are the experimental widths', bexp_widths, '\n #####################################################################')
    print('this was deemed the optimal alpha', optimal_alpha)
    print('these are the optimal alpha from least to most', optimal_alphas)

    plt.plot(alpha_arr, chi_squared_arr, color = 'royalblue')
    plt.ylabel(r"$\chi^2$")
    plt.xlabel(r"$\alpha$")
    plt.savefig('figs/chi_squared_bottonium.svg', bbox_inches = 'tight')
    plt.show()



'''bexp_widths = [ bex.onegtwog_percent*bex.total_width, bex.positron_percent*bex.total_width, bex.threegluons_percent*bex.total_width]
bexp_widths_errors = [ bex.onegtwog_error*bex.total_width, bex.positron_error*bex.total_width,bex.threegluons_error*bex.total_width]
btheo_widths_original = [bex.one_g_two_g(), bex.lepton_decay(), bex.three_gluons()]'''


def cex_chi_squared():

    cexp_widths = [cex.lepton_positron_percent *cex.total_width, cex.three_gluons_percent * cex.total_width, cex.one_g_two_g_percent * cex.total_width]
    ctheo_widths_original = [cex.lepton_decay(c.alpha_c), cex.three_gluons(c.alpha_c), cex.one_g_two_g(c.alpha_c)]

    #calculate chi_squared for every given chi
    type_c = 'CHARM'
    alpha_arr = values(type_c)
    chi_squared_arr = []
    for alpha in alpha_arr:
        print('the next alpha is being computed', alpha)
        #print('this is the first term of cex.lepton_decay', cex.lepton_decay(alpha))
        ctheo_widths = [cex.lepton_decay(alpha), cex.three_gluons(alpha), cex.one_g_two_g(alpha)]  #widths calculated for a given alpha
        if len(ctheo_widths) != len(cexp_widths):
            print('There are not as many experimental as theo widths! need to update! ')
        #sum_elements = [(theo-exp)**2/error**2 for theo, exp, error in zip(ctheo_widths, bexp_widths, bexp_widths_errors)] #
        #print('this is ctheo_widths', ctheo_widths)
        sum_elements = [(theo-exp)**2 for theo, exp in zip(ctheo_widths, cexp_widths)]
        chi_squared_value = sum(sum_elements)
        chi_squared_arr.append(chi_squared_value)

    idx_order = np.argsort(chi_squared_arr)
    arg_min_chi = idx_order[-1]
    #optimal_alpha = alpha_arr[arg_min_chi]
    optimal_alphas = [alpha_arr[i] for i in idx_order]
    optimal_alpha = optimal_alphas[0]

    ctheo_widths_optimal = [one_g_two_g(optimal_alpha),lepton_decay(optimal_alpha),three_gluons(optimal_alpha)]
    print('###################################################################### \n','These are the optimal theo widths', ctheo_widths_optimal, '\n #####################################################################')
    print('###################################################################### \n','These are the original theo widths', ctheo_widths_original, '\n #####################################################################')
    print('###################################################################### \n','These are the experimental widths', cexp_widths, '\n #####################################################################')
    print('this was deemed the optimal alpha', optimal_alpha)
    print('these are the optimal alpha from least to most', optimal_alphas)

    plt.plot(alpha_arr, chi_squared_arr, color = 'chocolate')
    plt.ylabel(r"$\chi^2$")
    plt.xlabel(r"$\alpha$")
    plt.savefig('figs/chi_squared_charmonium.svg', bbox_inches = 'tight')
    plt.show()






def c_ratio_chi_squared():
    #calculate chi_squared for every given chi
    c_ratio_exp_widths = [cex.total_width * cex.three_gluons_percent, cex.total_width * cex.one_g_two_g_percent]
    c_ratio_theo_widths_original = cex.ratios(c.alpha_c)

    type_c = 'CHARM'
    alpha_arr = values(type_c)
    chi_squared_arr = []
    for alpha in alpha_arr:
        print('the next alpha is being computed', alpha)
        #print('this is the first term of cex.lepton_decay', cex.lepton_decay(alpha))
        c_ratio_theo_widths = cex.ratios(alpha)  #widths calculated for a given alpha
        if len(c_ratio_theo_widths) != len(c_ratio_exp_widths):
            print('There are not as many experimental as theo widths! need to update! ')
        #sum_elements = [(theo-exp)**2/error**2 for theo, exp, error in zip(ctheo_widths, bexp_widths, bexp_widths_errors)] #
        #print('this is ctheo_widths', ctheo_widths)
        sum_elements = [(theo-exp)**2 for theo, exp in zip(c_ratio_theo_widths, c_ratio_exp_widths)]
        chi_squared_value = sum(sum_elements)
        chi_squared_arr.append(chi_squared_value)
    

    idx_order = np.argsort(chi_squared_arr)
    #optimal_alpha = alpha_arr[arg_min_chi]
    optimal_alphas = [alpha_arr[i] for i in idx_order]
    optimal_alpha = optimal_alphas[0]

    c_ratio_theo_widths_optimal = cex.ratios(optimal_alpha)
    print('###################################################################### \n','These are the optimal theo widths', c_ratio_theo_widths_optimal, '\n #####################################################################')
    print('###################################################################### \n','These are the original theo widths', c_ratio_theo_widths_original, '\n #####################################################################')
    print('###################################################################### \n','These are the experimental widths', c_ratio_exp_widths, '\n #####################################################################')
    print('these are the optimal alpha from least to most', optimal_alphas)

    plt.plot(alpha_arr, chi_squared_arr, color = 'chocolate')
    plt.ylabel(r"$\chi^2$")
    plt.xlabel(r"$\alpha$")
    plt.savefig('figs/ratios_chi_squared_charmonium.svg', bbox_inches = 'tight')
    plt.show()



def b_ratio_chi_squared():
    #calculate chi_squared for every given chi
    b_ratio_exp_widths = [bex.onegtwog_percent/bex.threegluons_percent, bex.threegluons_percent/bex.muon_percent]
    b_ratio_theo_widths_original = bex.ratios(c.alpha_b)

    type_b = 'BOTTOM'
    alpha_arr = values(type_b)
    chi_squared_arr = []
    for alpha in alpha_arr:
        print('the next alpha is being computed', alpha)
        #print('this is the first term of cex.lepton_decay', cex.lepton_decay(alpha))
        b_ratio_theo_widths = bex.ratios(alpha)  #widths calculated for a given alpha
        if len(b_ratio_theo_widths) != len(b_ratio_exp_widths):
            print('There are not as many experimental as theo widths! need to update! ')
        #sum_elements = [(theo-exp)**2/error**2 for theo, exp, error in zip(ctheo_widths, bexp_widths, bexp_widths_errors)] #
        #print('this is ctheo_widths', ctheo_widths)
        sum_elements = [(theo-exp)**2 for theo, exp in zip(b_ratio_theo_widths, b_ratio_exp_widths)]
        chi_squared_value = sum(sum_elements)
        chi_squared_arr.append(chi_squared_value)
    

    idx_order = np.argsort(chi_squared_arr)
    #optimal_alpha = alpha_arr[arg_min_chi]
    optimal_alphas = [alpha_arr[i] for i in idx_order]
    optimal_alpha = optimal_alphas[0]

    b_ratio_theo_widths_optimal = bex.ratios(optimal_alpha)
    print('###################################################################### \n','These are the optimal theo widths', b_ratio_theo_widths_optimal, '\n #####################################################################')
    print('###################################################################### \n','These are the original theo widths', b_ratio_theo_widths_original, '\n #####################################################################')
    print('###################################################################### \n','These are the experimental widths', b_ratio_exp_widths, '\n #####################################################################')
    print('this was deemed the optimal alpha', optimal_alpha)
    print('these are the optimal alpha from least to most', optimal_alphas)

    plt.plot(alpha_arr, chi_squared_arr, color = 'royalblue')
    plt.ylabel(r"$\chi^2$")
    plt.xlabel(r"$\alpha$")
    plt.savefig('figs/ratios_chi_squared_bottonium.svg', bbox_inches = 'tight')
    plt.show()


cex_chi_squared()
#bex_chi_squared()
#c_ratio_chi_squared()
#b_ratio_chi_squared()
#make it such that it prints optimal values for the best alpha/chi-squared



#optimal_alpha_b = 0.17131313131313133
#optimal_alpha_c = 0.26121212121212123 