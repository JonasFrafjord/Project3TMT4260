import sys
import numpy as np
import math
from matplotlib import pyplot as plt

#Readme:
    #subscript t, i.e. _t, means temporary, of functional variable. It lives inside the scope of the function (most times at least)



#Global variables:
T_k = 273.15            #Celcius to Kelvin, 1 celcius
T_m = 660.5#+T_k         #Melt temp pure Al, Kelvin
T_e = 577.0#+T_k         #Eutectic temp Al-Si, Kelvin
C_s = 1.5               #Solubility of Si at T_e, wt%Si
C_e = 12.2              #Eutectic composition, wt%Si
k_pc = C_s/C_e          #Partitioning coefficient defined to be C_sol/C_liq, is constant due to linearised phase diagram
m_upper = (T_m-T_e)/C_e #rate of linear line sol-liq
m_lower = (T_m-T_e)/C_s #rate of linear line sol-sol

f1, subfig1 = plt.subplots(1,2, sharey=True)
plt.suptitle('Weight fraction of solid as a function of temperature')
plt.ylim(-0,1.01)
subfig1[0].set_title('1wt%')
subfig1[1].set_title('8wt%')
subfig1[0].set_ylabel('f')
subfig1[0].set_xlabel('T[C]')
subfig1[1].set_xlabel('T[C]')
f2, subfig2 = plt.subplots(1,2, sharey=True)
subfig2[0].set_title('1wt%')
subfig2[1].set_title('8wt%')
subfig2[0].set_ylabel('df/dt')
subfig2[0].set_xlabel('T[C]')
subfig2[1].set_xlabel('T[C]')

plt.suptitle('Weight fraction of solid differentiated with respect to temperature')

#The temperature associated to a given concentration C_0, not a free variable-sol-liq
def getT_L(C_0_t):
    return T_m-m_upper*C_0_t
#The temperature associated to a given concentration C_liq_t, a free variable of the system
def getT(C_liq_t):
    return T_m-m_upper*C_liq_t
#The temperature associated to a given concentration C_0, not a free variable, sol-sol
def getT_S(C_0_t):
    return T_m-m_lower*C_0_t

#The weight fraction of solid at equilibrium, given by the lever rule
def SF_Equi(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_t<T_S_t: return 1 #All is solid, do not distinguish between diferent solid phases
    return 1/(1-k_pc)*(T_L_t-T_t)/(T_m-T_t)
#The weight fraction differentiated with respect to temperature
def SF_Equi_dt(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_t<T_S_t: return 0 #All is solid, do not distinguish between diferent solid phases
    return 1/(k_pc-1)*(T_m-T_L_t)/(T_m-T_t)**2

#The weight fraction of solid in the Scheil model
def SF_scheil(T_L_t, T_S_t, T_t):
#    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
#    if T_t<T_S_t: return 1 #All is solid, do not distinguish between diferent solid phases
    if T_m == T_t:return 0
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))> 1:return 1
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))< 0:return 0
    return 1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1))
#The weight fraction differentiated with respect to temperature
def SF_scheil_dt(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    #if T_t<T_S_t: return 0 #All is solid, do not distinguish between diferent solid phases
    if T_m == T_L_t: return 1
    return (1/(k_pc-1))*(1/(T_m-T_L_t))*((T_m-T_t)/(T_m-T_L_t))**((2-k_pc)/(k_pc-1))


def homework2():
    Nt = 1e3
    T_min = T_e
    T_max = T_m
    C_0 = [1.0, 8.0]
    T_L = [getT_L(i) for i in C_0]
    T_S = [getT_S(i) for i in C_0]
    T = np.linspace(T_min, T_max, Nt)
    F_s_eq = [[SF_Equi(T_L[0],T_S[0], j) for j in T],[SF_Equi(T_L[1],T_S[1], j) for j in T]]
    F_s_eq_dt = [[SF_Equi_dt(T_L[0],T_S[0], j) for j in T],[SF_Equi_dt(T_L[1],T_S[1], j) for j in T]]
    subfig1[0].plot(T,F_s_eq[0], label = 'Equilibrium')
    subfig1[1].plot(T,F_s_eq[1], label = 'Equilibrium')
    subfig2[0].plot(T,F_s_eq_dt[0], label = 'Equilibrium')
    subfig2[1].plot(T,F_s_eq_dt[1], label = 'Equilibrium')
    print('The liquid fraction at the eutectic temperature is 0 for {} wt%Si, assuming equilibrium.'.format(C_0[0]))
    print('The liquid fraction at the eutectic temperature is {0:.3f} for {1} wt%Si, assuming equilibrium.'.format(1-SF_Equi(T_L[1],T_S[1],T_e),C_0[1]))

def homework3():
    Nt = 1e3
    T_min = T_e
    T_max = T_m
    C_0 = [1.0, 8.0]
    T_L = [getT_L(i) for i in C_0]
    T_S = [getT_S(i) for i in C_0]
    T = np.linspace(T_min, T_max, Nt)
    F_s_sch = [[SF_scheil(T_L[0],T_S[0], j) for j in T],[SF_scheil(T_L[1],T_S[1], j) for j in T]]
    F_s_sch_dt = [[SF_scheil_dt(T_L[0],T_S[0], j) for j in T],[SF_scheil_dt(T_L[1],T_S[1], j) for j in T]]
    subfig1[0].plot(T,F_s_sch[0], ls='--', label = 'Scheil')
    subfig1[1].plot(T,F_s_sch[1], ls='--', label = 'Scheil')
    subfig2[0].plot(T,F_s_sch_dt[0], ls='--', label = 'Scheil')
    subfig2[1].plot(T,F_s_sch_dt[1], ls='--', label = 'Scheil')
    print('The liquid fraction at the eutectic temperature is {0:.3f} for {1} wt%Si in the Scheil model.'.format(1-SF_scheil(T_L[0],T_S[0],T_e),C_0[0]))
    print('The liquid fraction at the eutectic temperature is {0:.3f} for {1} wt%Si in the Scheil model.'.format(1-SF_scheil(T_L[1],T_S[1],T_e),C_0[1]))

def main(argv):
    homework2()
    homework3()
    subfig1[1].legend()
    subfig2[1].legend(loc='best')
    plt.show()

#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
