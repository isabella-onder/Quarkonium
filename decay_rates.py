import scipy as sp
import numpy as np
import matplotlib.pyplot as plt
import constants as c
import cextension as cex
import bextension as bex
 
#just printing all the results: both numerical and experimental
alpha_c = c.alpha_c
alpha_c = 0.26121212121212123 
all_widths_charmonium = [cex.lepton_decay(alpha_c), cex.three_gluons(alpha_c), cex.one_g_two_g(alpha_c), cex.alternative_mag_transition(1), cex.three_photons()]
all_ratios_charmonium = cex.ratios(c.alpha_c)
experimental_widths_charmonium = [cex.lepton_positron_percent *cex.total_width, cex.three_gluons_percent * cex.total_width, cex.one_g_two_g_percent * cex.total_width, cex.hyperfine_percent * cex.total_width, cex. three_photons_percent * cex.total_width]
experimental_ratios_charmonium = [cex.total_width * cex.three_gluons_percent, cex.total_width * cex.one_g_two_g_percent]

'''
print('#################################################\n these are all my numerical widths \n lepton, three gluon, one photon two gluons, spin flip, three photons')
print(all_widths_charmonium, '\n#####################################################')

print('#################################################\nthese are the corresponding experimental widths \n lepton, three gluon, one photon two gluons, spin flip. three photons')
print(experimental_widths_charmonium, '\n ##############################################################################')

print('and for the ratios of three gluons and one gluon two photons', all_ratios_charmonium, '\n', experimental_ratios_charmonium)
'''

alpha_b = c.alpha_b
alpha_b = 0.17131313131313133
experimental_widths_bottonium = [ bex.onegtwog_percent*bex.total_width, bex.positron_percent*bex.total_width, bex.threegluons_percent*bex.total_width]
experimental_ratios_bottonium = [bex.onegtwog_percent/bex.threegluons_percent, bex.threegluons_percent* bex.muon_percent]
all_widths_bottonium = [bex.one_g_two_g(alpha_b), bex.lepton_decay(alpha_b), bex.three_gluons(alpha_b)]
all_ratios_bottonium = bex.ratios(c.alpha_b)
bottonium_text_widths = 'one photon two gluons, lepton, three gluons'
bottonium_text_ratios = 'one photon two gluons, muons'

charmonium_text_widths = 'lepton, three gluon, one photon two gluons, spin flip, three photons'
charmonium_text_ratios = 'ratios of three gluons and one gluon two photons'
#do fractional differences to get an idea of how right or wrong
def fractional_difference(numerical_widths, experimental_widths, text_order):
    fractional_differences = []
    for num_width, exp_width,  in zip(numerical_widths, experimental_widths):
        frac_diff = np.abs((num_width - exp_width)/exp_width)
        fractional_differences.append(frac_diff)
    print('these are the fractional differences of \n', text_order, '\n', fractional_differences)
    return fractional_differences 

#fractional_difference(all_widths_charmonium, experimental_widths_charmonium, charmonium_text_widths)
#print(all_widths_charmonium)
#print(experimental_widths_charmonium)
#fractional_difference(all_ratios_charmonium, experimental_ratios_charmonium, charmonium_text_ratios)

#fractional_difference(all_widths_bottonium, experimental_widths_bottonium, bottonium_text_widths)
print(all_widths_bottonium)
print(experimental_widths_bottonium)
#fractional_difference(all_ratios_bottonium, experimental_ratios_bottonium, bottonium_text_ratios)
