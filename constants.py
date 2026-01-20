#assume alpha is fixed and beta is already found from the beta file
alpha_c = 0.40
alpha_b = 0.28

beta = 0.210083

m_c = 1.2730 #Gev.c^-2
m_b = 4.183

rmax = 30 #because it has been sufficient for now




#unit converters:
def kg_to_GeVc2(m):                 #input mass in kg
    return m * 5.610 * 10 **26


def m_to_GeV(r):                    #input length in m
    return r * 5.0677 * 10 ** 15