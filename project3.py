#!/usr/bin/env python3

#---------------------------------------
#Team 1:
#Tina Bergh
#Jonas Frafjord
#Jonas K. Sunde
#---------------------------------------

import sys
import numpy as np
import math
from matplotlib import pyplot as plt

#Readme:
    #subscript t, i.e. _t, means temporary, of functional variable. It lives inside the scope of the function (most times at least)
    #The main function will run when you execute the script. Here you can change which functions will run, and which parameters you want to set

#Required information:
	#Nucleation temperature is assumed to be 0.1 degrees below the melting temperature, i.e. T_L-0.1
	#The referance temperature is chosen to be 2 degrees below the melting temperature, i.e. T_L-2

#Global variables:
T_k = 273.15            #Celcius to Kelvin, 1 celcius
T_m = 660.5+T_k         #Melt temp pure Al, Kelvin
T_e = 577.0+T_k         #Eutectic temp Al-Si, Kelvin
C_s = 1.5               #Solubility of Si at T_e, wt%Si
C_e = 12.2              #Eutectic composition, wt%Si
C_0r = 4				#Chosen referance concentration in wt%Si, starting concentration
X_c						#The MoleFraction it took t_s (i.e. t_star) seconds to form
k_pc = C_s/C_e          #Partitioning coefficient defined to be C_sol/C_liq, is constant due to linearised phase diagram
m_upper = (T_m-T_e)/C_e #rate of linear line sol-liq
m_lower = (T_m-T_e)/C_s #rate of linear line sol-sol
L = 0.8					#Latent heat J/mm^3
rhoC = 0.0027			#Volume heat capacity [J/(Celcius*mm^3)] density dependent
lambdaL = 0.095			#Thermal conductivity liquid [W/(Celcius*mm)]
lambdaS = 0.21			#Thermal conductivity solid [W/(Celcius*mm)]
t_r = 6					#6 seconds simulation

	
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
#The referance temperature
def getT_r(T_L_t):
	return T_L_t - 2

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

#The solid weight fraction using the Scheil model. Non-equilibrium freezing.
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
    



def main(argv):
    print('Program started')

#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
