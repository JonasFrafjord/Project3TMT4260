#!/usr/bin/env python3

#---------------------------------------
#Group 1:
#Tina Bergh
#Jonas Frafjord
#Jonas K. Sunde
#---------------------------------------

import sys
import numpy as np
import math
from matplotlib import pyplot as plt

#Readme:
    #Subscript t, i.e. _t, means temporary, of functional variable. It lives inside the scope of the function (most times at least)
    #The main function will run when you execute the script. Here you can change which functions will run, and which parameters you want to set

#Required information:
	#Nucleation temperature is assumed to be 0.1 degrees below the melting temperature, i.e. T_L-0.1
	#The reference temperature is chosen to be 2 degrees below the melting temperature, i.e. T_L-2.0

#Global variables:
T_k = 273.15            #Degrees Kelvin at 0 degrees Centigrade
T_m = 660.5#+T_k         #Melting temperature of pure Al [K]
T_e = 577.0#+T_k         #Eutectic temperature of binary Al-Si [K]
C_s = 1.5               #Solubility of Si in Al at T_e, wt% Si
C_e = 12.2              #Eutectic composition, wt% Si

#C_0r = 				#Chosen reference concentration in wt% Si, starting concentration
#X_c = 4.0					#The MoleFraction it took t_s (i.e. t_star) seconds to form

k_pc = C_s/C_e          #Partitioning coefficient defined to be C_sol/C_liq, is constant due to linearised phase diagram
m_upper = (T_m-T_e)/C_e #rate of linear line sol-liq
m_lower = (T_m-T_e)/C_s #rate of linear line sol-sol
#L = 0.8					#Latent heat [J/mm^3]
L = 1.0746				#Latent heat [J/mm^3]
rho = 2.7e-3           #Density[g/mm^3]
C_hc = 24.20            #Heat capacity [J/mol]
M_mAl = 26.98          #[g/mol] Molar mass Aluminium
rhoC = rho/M_mAl*C_hc #Volume heat capacity [J/(Celcius*mm^3)]
print('rhoC1 = {0:.5f} J/Cmm^3\n'.format(rhoC2))
rhoC2 = 0.0027			   #Volume heat capacity [J/(Celcius*mm^3)] density dependent
print('rhoC2 = {0:.5f} J/Cmm^3'.format(rhoC))
#lambdaL = 0.095			#Thermal conductivity liquid Al [W/(Celcius*mm)]
#lambdaS = 0.21			#Thermal conductivity solid Al [W/(Celcius*mm)]
#Andre verdiar frå literaturen:
lambdaL = 0.094			#Thermal conductivity liquid Al [W/(Celcius*mm)] @ 665 degrees Celcius
lambdaS = 0.213			#Thermal conductivity solid Al [W/(Celcius*mm)] @ 630 degrees Celcius
t_r = 6.0					#6 seconds simulation


	
#The temperature associated to a given concentration C_0. Not a free variable, but given by sol-liq-line.
def getT_L(C_0_t):
    return T_m-m_upper*C_0_t
#The temperature associated to a given concentration C_liq_t, a free variable of the system
def getT(C_liq_t):
    return T_m-m_upper*C_liq_t
#The temperature associated to a given concentration C_0. Not a free variable, but given by sol-sol-line.
def getT_S(C_0_t):
    return T_m-m_lower*C_0_t
#The nucleation temperature. We may delete this function and the next one, at some point. If we find it redundant.
def getT_n(T_L_t):
	return T_L_t - 0.1
 
#The solid weight fraction at equilibrium, given by the lever rule.
def SF_Equi(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_t<T_S_t: return 1 #All is solid, do not distinguish between diferent solid phases
    return 1/(1-k_pc)*(T_L_t-T_t)/(T_m-T_t)
#The solid weight fraction differentiated with respect to temperature
def SF_Equi_dt(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_t<T_S_t: return 0 #All is solid, do not distinguish between diferent solid phases
    return 1/(k_pc-1)*(T_m-T_L_t)/(T_m-T_t)**2

#The solid weight fraction using the Scheil model. Non-equilibrium solidification.
def SF_scheil(T_L_t, T_S_t, T_t):
    if T_m == T_t:return 0
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))> 1:return 1
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))< 0:return 0
    return 1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1))
#The solid weight fraction differentiated with respect to temperature
def SF_scheil_dt(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_m == T_L_t: return 1
    return (1/(k_pc-1))*(1/(T_m-T_L_t))*((T_m-T_t)/(T_m-T_L_t))**((2-k_pc)/(k_pc-1))
    
#Choice of reference state:
X_c = 0.15					#Reference MoleFraction. The MoleFraction it took t_s (i.e. t_star) seconds to form
C_0r = X_c*C_s#SF_Equi()	#Resulting reference concentration in wt% Si
#The reference temperature
def getT_r(T_L_t):
	return T_L_t - 2.0
	
#Mole fraction as a function of a referance mole fraction and time, X_c and t_s respectively, time t and the integer n.
def XMF(X_c_t, n, t_t, t_s_t=1.0):
    return 1-(1-X_c_t)**((t_t/t_s_t)**n)
#Differentiate numerically XMF as a function of time, by 1st order finite difference method.
def dXMFdt(X_c_plus_t, X_c_minus_t, dt_t):
    return (X_c_plus_t-X_c_minus_t)/(2*dt_t)
#An analytic solution for dXMF/dt
def dXMFdt_anal(X_t, X_c_t, n, t_s_t=1.0):
    if X_t==0: return 0
    return -(1-X_t)*math.log(1-X_t)*n/t_s_t/(math.log(1-X_t)/math.log(1-X_c_t))**(1/n)

#Precursor to differential form of MoleFraction
#Kick-off equation
def dXMFdt_precursor(X_c_t, n, t_t,t_s_t):
    if X_c_t==0: return 0
    return -(1-X_c_t)**((t_t/t_s_t)**n)*math.log(1-X_c_t)*((t_t/t_s_t)**n)*n/t_t #Solve as backward Euler

#Time constant
#NB: Define all parameters (calc from chosen reference condition)
def t_star(Delta_T_temp,C_0_temp,N_temp,f_s_temp,n):
    return t_ref*(Delta_T_ref/Delta_T_temp)**2*(C_0_temp/C_0_ref)*(N_r/N_temp)**(1/n)*(f_s_temp/f_s_r)**(1/n)

def solidification():
    #Temperature boundaries
    T_min = T_e
    T_max = T_m
    
    #Temporal discretisation
    Nt = int(1e3) 
    tmax = 6.0 #Seconds of simulation
    dt = tmax/Nt #Time increments
    t = np.linspace(0,tmax,Nt)
    X_c = [0.05, 0.15] #Reference MoleFraction
    
    #Kick-off MoleFraction
    #X_0 = dXMFdt_precursor(...)
    
    #Equilibrium solid fraction
    C_0 = [1.0, 8.0] # wt% Si
    T_L = [getT_L(i) for i in C_0]
    [print('T_L = {0:.2f} deg C for C_0 = {1:.2f} wt% Si'.format(getT_L(i),i)) for i in C_0]
    T_S = [getT_S(i) for i in C_0]
    T = np.linspace(T_min, T_max, Nt)
    F_s_eq = [[SF_Equi(T_L[0],T_S[0], j) for j in T],[SF_Equi(T_L[1],T_S[1], j) for j in T]]
    F_s_eq_dt = [[SF_Equi_dt(T_L[0],T_S[0], j) for j in T],[SF_Equi_dt(T_L[1],T_S[1], j) for j in T]]
    
    #Scheil model 'non-eq. lever rule'
    F_s_sch = [[SF_scheil(T_L[0],T_S[0], j) for j in T],[SF_scheil(T_L[1],T_S[1], j) for j in T]]
    F_s_sch_dt = [[SF_scheil_dt(T_L[0],T_S[0], j) for j in T],[SF_scheil_dt(T_L[1],T_S[1], j) for j in T]]    
    
    #MoleFraction evolution
    n = [1,2,3]
    nlist = np.append(n,n)
    X_c_n = [[X_c[0],j] for j in n] + [[X_c[1],j] for j in n]
    X = [[XMF(i,j,k) for k in t] for i in X_c for j in n]
    dXdt_anal = [[dXMFdt_anal(i[k], j[0], j[1]) for k in range(Nt)] for i,j in zip(X,X_c_n)]
    #plt.figure()
    dXdt = [[dXMFdt(Xlist[k+1],Xlist[k-1],tmax/Nt) for k in range(1,Nt-1)] for Xlist in X]

    #Main loop
    for i in t:
        # Calc Kick-off value X_0 from dXMFdt_precursor(..)
        # --> Calc. X_next from dXMFdt_anal(...)        
        # --> Calc. f_s_next from df_s/dt = f_m*dX/dt
        # --> Calc T_next from dT/dt = -dotQ/rhoC+L/rhoC*df_s/dt
def main(argv):
    print('Program started')
    solidification()
    
#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
