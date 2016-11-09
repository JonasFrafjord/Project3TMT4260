import sys
import numpy as np
import math
from matplotlib import pyplot as plt

#Readme:
    #subscript t, i.e. _t, means temporary, of functional variable. It lives inside the scope of the function (most times at least)
    #The main function will run when you execute the script. Here you can change which functions will run


#Global variables:
T_k = 273.15            #Celcius to Kelvin, 1 celcius
T_m = 660.5#+T_k         #Melt temp pure Al, Kelvin
T_e = 577.0#+T_k         #Eutectic temp Al-Si, Kelvin
C_s = 1.5               #Solubility of Si at T_e, wt%Si
C_e = 12.2              #Eutectic composition, wt%Si
k_pc = C_s/C_e          #Partitioning coefficient defined to be C_sol/C_liq, is constant due to linearised phase diagram
m_upper = (T_m-T_e)/C_e #rate of linear line sol-liq
m_lower = (T_m-T_e)/C_s #rate of linear line sol-sol


#This is a setup for the figures which we will plot on. The plots are added when we need to.
#Will be a weight fraction plot
if False:
    f1, subfig1 = plt.subplots(1,2, sharey=True) 
    plt.suptitle('Weight fraction of solid as a function of temperature')
    plt.ylim(-0,1.01)
    subfig1[0].set_title('1wt%')
    subfig1[1].set_title('8wt%')
    subfig1[0].set_ylabel('f')
    subfig1[0].set_xlabel('T[C]')
    subfig1[1].set_xlabel('T[C]')
# Will contain the differentiated weight fraction
if False:
    f2, subfig2 = plt.subplots(1,2, sharey=True)
    subfig2[0].set_title('1wt%')
    subfig2[1].set_title('8wt%')
    subfig2[0].set_ylabel('df/dt')
    subfig2[0].set_xlabel('T[C]')
    subfig2[1].set_xlabel('T[C]')
    plt.suptitle('Weight fraction of solid differentiated with respect to temperature')
# HW4. Mole fraction
if False:
    f3, subfig3 = plt.subplots(1,2)
    subfig3[0].set_title('1wt%')
    subfig3[1].set_title('8wt%')
    subfig3[0].set_ylabel('df/dt')
    subfig3[0].set_xlabel('T[C]')
    subfig3[1].set_xlabel('T[C]')
    plt.suptitle('Weight fraction of solid differentiated with respect to temperature')
if True:
    f4, subfig4 = plt.subplots(1,2)
    subfig4[0].set_title('1wt%')
    subfig4[1].set_title('8wt%')
    subfig4[0].set_ylabel('df/dt')
    subfig4[0].set_xlabel('T[C]')
    subfig4[1].set_xlabel('T[C]')
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
#Mole fraction as a function of a referance mole fraction and time (X_c, t_s) respectively, time t and the integer n.
def XMF(X_c_t, n, t_t, t_s_t=1.0):
    return 1-(1-X_c_t)**((t_t/t_s_t)**n)
#Differentiate numerically XMF as a function of time, next incriment
def dXMFdt(X_c_plus_t, X_c_minus_t, dt_t):
    return (X_c_plus_t-X_c_minus_t)/(2*dt_t)
def dXMFdt_anal(X_t, X_c_t, n, t_s_t=1.0):
    if X_t==0: return 0
    return -(1-X_t)*math.log(1-X_t)*n/t_s_t/(math.log(1-X_t)/math.log(1-X_c_t))**(1/n)
    


#Equilibrium solution
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

#Scheil model
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

def homework4():
    Nt = int(1e3)
    t = np.linspace(0,5,Nt)
    X_c = [0.05, 0.15]
    n = [1,2,3]
    X = [[XMF(i,j,k) for k in t] for i in X_c for j in n]
    #plt.figure()
    dXdt = [[dXMFdt(Xlist[k+1],Xlist[k-1],1/Nt) for k in range(1,Nt-1)] for Xlist in X]
    nlist = np.append(n,n)
    for Xlist,n_id in zip(X,nlist):
        subfig3[0].plot(t,Xlist, label='{}'.format(n_id))
    #plt.figure()
    for Xlist, n_id in zip(dXdt,nlist):
        subfig3[1].plot(t[1:Nt-1],Xlist, label='{}'.format(n_id))

def homework5():
    Nt = int(1e3)
    t = np.linspace(0,5,Nt)
    X_c = [0.05, 0.15]
    n = [1,2,3]
    nlist = np.append(n,n)
    X_c_n = [[X_c[0],j] for j in n] + [[X_c[1],j] for j in n]
    X = [[XMF(i,j,k) for k in t] for i in X_c for j in n]
    dXdt_anal = [[dXMFdt_anal(i[k], j[0], j[1]) for k in range(Nt)] for i,j in zip(X,X_c_n)]
    for Xlist,n_id in zip(X,nlist):
        subfig4[0].plot(t,Xlist, label='{}'.format(n_id))
    for Xlist, n_id in zip(dXdt_anal,nlist):
        subfig4[1].plot(t,Xlist, label='{}'.format(n_id))

def main(argv):
    #homework2()
    #homework3()
#    homework4()
    homework5()
#    subfig1[1].legend()
#    subfig2[1].legend(loc='best')
#    subfig3[1].legend()
    subfig4[1].legend()
    plt.show()

#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
