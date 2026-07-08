#assume alpha is fixed and beta is already found from the beta file


#for when trying to do it with the actual strong coupling constant value
#for the energy range of our charmonia
#beta = 0.143457
#alpha_c = 0.225
#for trying to do with just two charms i.e. proper lower range and as measured
#alpha_c = 0.26
#beta = 0.15458984
#the ones I had been using incorrectly for a while
#alpha_c = 0.40
#beta = 0.210083



#with the - error
#alpha_c = 0.35
#beta = 0.18809570312499996

#with the + error
#alpha_c = 0.41
#beta = 0.21475585937500002

#for the literature value .38
alpha_c = 0.38
beta_c = 0.20087890625000002 #kept the approx with mc = 1.27

alpha_b = 0.28
beta_b = 0.209970703125

alpha_b = 0.33
beta_b = 0.21940429687500002



m_c = 1.2730 #Gev.c^-2 p/m 0.0028
m_b = 4.183

rmax = 30 #because it has been sufficient for now

e_c = 2/3 #charm charge
e_b = -1/3 #bottom charge


#unit converters:
def kg_to_GeVc2(m):                 #input mass in kg
    return m * 5.610 * 10 **26


def m_to_GeV(r):                    #input length in m
    return r * 5.0677 * 10 ** 15