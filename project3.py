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
print('rhoC1 = {0:.5f} J/Cmm^3'.format(rhoC2))
rhoC2 = 0.0027			   #Volume heat capacity [J/(Celcius*mm^3)] density dependent
print('rhoC2 = {0:.5f} J/Cmm^3'.format(rhoC))
#lambdaL = 0.095			#Thermal conductivity liquid Al [W/(Celcius*mm)]
#lambdaS = 0.21			#Thermal conductivity solid Al [W/(Celcius*mm)]
#Andre verdiar frå literaturen:
lambdaL = 0.094			#Thermal conductivity liquid Al [W/(Celcius*mm)] @ 665 degrees Celcius
lambdaS = 0.213			#Thermal conductivity solid Al [W/(Celcius*mm)] @ 630 degrees Celcius
t_sim = 6.0					#6 seconds simulation

#This is a setup for the figures which we will plot on. The plots are added when we need to.
#Will be a weight fraction plot.
if False:
    f1, subfig1 = plt.subplots(1,2)
    #f3.tight_layout()
    subfig1[0].set_title(r'X$_{c}$ = 0.15',y=1.01)
    subfig1[1].set_title(r'X$_{c}$ = 0.15',y=1.01)
    subfig1[0].set_ylabel(r'f$_{s}$')
    subfig1[0].set_xlabel('t [s]')
    subfig1[1].set_xlabel('t [s]')
    subfig1[1].set_ylabel(r'df$_{s}$/dt')
    plt.suptitle(r'Weight fraction of solid $\alpha$-phase (left) and time derivative (right)', position=(0.7,1.05),fontsize= 18)
    
    left  = 0.125  # the left side of the subplots of the figure
    right = 1.2    # the right side of the subplots of the figure
    bottom = 0.1   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.3   # the amount of width reserved for blank space between subplots
    hspace = 0.5   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
    
if False:
    f2, subfig2 = plt.subplots(4,1)
    #f3.tight_layout()
    subfig2[0].set_title(r'Temp',y=1.01)
    subfig2[0].set_ylabel(r'T [$^{\circ}$C]')
    subfig2[0].set_xlabel('t [s]')
    subfig2[1].set_title(r'Scaled volume fraction',y=1.01)
    subfig2[1].set_ylabel(r'X')
    subfig2[1].set_xlabel('t [s]')
    subfig2[2].set_title(r'Volume fraction $\alpha$',y=1.01)
    subfig2[2].set_ylabel(r'f$_{s}$')
    subfig2[2].set_xlabel('t [s]')
    subfig2[3].set_title(r'Time constant t_star',y=1.01)
    subfig2[3].set_ylabel(r't$_{star}$')
    subfig2[3].set_xlabel('#')
    plt.suptitle(r'Solidification of $\alpha$-phase (test-plot)', position=(0.7,1.05),fontsize= 18)
    
    left  = 0.125  # the left side of the subplots of the figure
    right = 1.2    # the right side of the subplots of the figure
    bottom = 0.3   # the bottom of the subplots of the figure
    top = 0.9      # the top of the subplots of the figure
    wspace = 0.3   # the amount of width reserved for blank space between subplots
    hspace = 0.5   # the amount of height reserved for white space between subplots
    plt.subplots_adjust(left, bottom, right, top, wspace, hspace)
	
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

###################################################
#Choice of reference state:
#Defined from t = t_r, T=T_r and N = N_r
t_r = 1.0                  #Reference time, arb. chosen
N_r = 1000                  #Reference nucleation sites, arb. chosen
#The reference temperature
def getT_r(T_L_t):
	return T_L_t - 2.0

# X = X_c when t = t_star
X_c = 0.15					#Reference MoleFraction. The MoleFraction it took t_s (i.e. t_star) seconds to form
C_0r = X_c*C_s#SF_Equi()	#Resulting reference concentration in wt% Si
f_s_r = 2.0                #Reference volume fraction of spherical precipitates

#Time constant t_star (Applying the Hunt model, i.e. growth rate V prop. to undercooling^2/C_0
#NB: Define all parameters (calc from chosen reference condition)
def get_t_star(T_L_temp,T_temp,C_0_temp,N_temp,f_s_temp,n):
    return t_r*((getT_r(T_L_temp)-T_temp)/(T_L_temp-T_temp))**2*(C_0_temp/C_0r)*(N_r/N_temp)**(1/n)*(f_s_temp/f_s_r)**(1/n)

####################################################
	
#Mole fraction as a function of a referance mole fraction and time, X_c and t_s respectively, time t and the integer n.
def XMF(X_c_t, n, t_t, t_s_t=1.0):
    return 1-(1-X_c_t)**((t_t/t_s_t)**n)
#Differentiate numerically XMF as a function of time, by 1st order finite difference method.
def dXMFdt(X_c_plus_t, X_c_minus_t, dt_t):
    return (X_c_plus_t-X_c_minus_t)/(2*dt_t)
#An analytic solution for dXMF/dt
def dXMFdt_anal(X_temp,X_c_temp,n,t_temp,t_star_temp):
    if X_temp==0: return 0
    return -(1-X_temp)*math.log(1-X_temp)*n/t_star_temp/((math.log(1-X_temp)/math.log(1-X_c_temp))**(1/n))

#Precursor to differential form of MoleFraction
#Kick-off equation
def dXMFdt_precursor(X_c_temp, n, t_temp,t_star_temp,dt_temp,X_start_temp):
    if X_c_temp==0: return 0
    return X_start_temp-dt_temp*(1-X_c_temp)**((t_temp/t_star_temp)**n)*math.log(1-X_c_temp)*((t_temp/t_star_temp)**n)*n/t_temp #Solve as backward Euler

#Total solidification time, eq. 15
def sol_time(t_eq_temp,f_eq_temp,L_temp,a_max_temp,rhoC_temp):
    return t_eq_temp-f_eq_temp*L_temp/a_max_temp/rhoC_temp
    
#Intrinsic cooling rate during solidification
def dT_next(dotQ_temp,rhoC_temp,L_temp,dfdT_Scheil_temp,dt_temp,T_prev_temp):
    return T_prev_temp-dt_temp*dotQ_temp/rhoC_temp+dt_temp*L_temp/rhoC_temp*dfdT_Scheil_temp

#Steady state cooling rate, i.e. dynamic balance occuring when X --> 1
def dT_next_steady_state(L_temp,rhoC_temp,dfdT_Scheil,a_max_temp):
    return a_max_temp*(L_temp/rhoC_temp*dfdT_Scheil-1)**-1

def solidification():
    #Temperature boundaries
    T_min = T_e
    T_max = T_m
    #Nucleation sites
    N = 1000
    C_0 = [1.0, 4.0] # wt% Si, Alloy composition
    
    #Temporal discretisation
    Nt = int(1e3) 
    tmax = 6.0 #Seconds of simulation
    dt = tmax/Nt #Time increments
    t = np.linspace(0,tmax,Nt)
    X_c = [0.05, 0.15] #Reference MoleFraction
    
    #Kick-off MoleFraction
    #X_0 = dXMFdt_precursor(...)
    
    #Equilibrium solid fraction
    T_L = [getT_L(i) for i in C_0]
    [print('T_L = {0:.2f} deg C for C_0 = {1:.2f} wt% Si'.format(getT_L(i),i)) for i in C_0]
    T_S = [getT_S(i) for i in C_0]
    
    T_test = [T_max-15,T_max-50]
    #Solid fraction
    F_s_eq = [SF_Equi(T_L[0],T_S[0],T_test[0]),SF_Equi(T_L[1],T_S[1], T_test[1])]
    print(F_s_eq)
    [print('F_s_eq{0:d} = {1:.3f} wt% Si'.format(i+1,F_s_eq[i])) for i in range(2)]
    
    #t_star
    t_star = [get_t_star(T_L[i],T_test[i],C_0[i],N,F_s_eq[i],n=3) for i in range(2)]
    [print('t_star{0:d} = {1:.3f} s'.format(i+1,t_star[i])) for i in range(2)]


    T = np.linspace(T_min, T_max, Nt)
    #F_s_eq = [[SF_Equi(T_L[0],T_S[0], j) for j in T],[SF_Equi(T_L[1],T_S[1], j) for j in T]]
    #F_s_eq_dt = [[SF_Equi_dt(T_L[0],T_S[0], j) for j in T],[SF_Equi_dt(T_L[1],T_S[1], j) for j in T]]
    
    #Scheil model 'non-eq. lever rule'
    #F_s_sch = [[SF_scheil(T_L[0],T_S[0], j) for j in T],[SF_scheil(T_L[1],T_S[1], j) for j in T]]
    #F_s_sch_dt = [[SF_scheil_dt(T_L[0],T_S[0], j) for j in T],[SF_scheil_dt(T_L[1],T_S[1], j) for j in T]]    
    
    #MoleFraction evolution
    #n = [1,2,3]
    #nlist = np.append(n,n)
    #X_c_n = [[X_c[0],j] for j in n] + [[X_c[1],j] for j in n]
    #X = [[XMF(i,j,k) for k in t] for i in X_c for j in n]
    #dXdt_anal = [[dXMFdt_anal(i[k], j[0], j[1]) for k in range(Nt)] for i,j in zip(X,X_c_n)]
    #plt.figure()
    #dXdt = [[dXMFdt(Xlist[k+1],Xlist[k-1],tmax/Nt) for k in range(1,Nt-1)] for Xlist in X]


    # Calc Kick-off value X_0 from dXMFdt_precursor(..)
    X_0 = [dXMFdt_precursor(X_c[i], 3, 1,t_star[i],dt,F_s_eq[i]) for i in range(2)]
    print(X_0)
    
    # Iteration steps
    sim = 10
    X = np.zeros(sim)
    X[0] = X_0[0]
    f_s = np.zeros(sim)
    f_s[0] = F_s_eq[0]
    T = np.zeros(sim)
    T[0] = T_test[0]
    t_st = np.zeros(sim)
    t_st[0] = t_star[0]
    n_=3
    
    #Main loop
    #for i in t:
    for i in range(sim-1):
        # --> Calc. X_next from dXMFdt_anal(...)
        t_curr = i*dt
        X[i+1] = dXMFdt_anal(X[i],X_c[0],n_,t_curr,t_st[i])
        # --> Calc. f_s_next from df_s/dt = f_m*dX/dt
        f_s[i+1] = C_s*X[i+1] # NB test with C_s, should be f_m
        # --> Calc T_next from dT/dt = -dotQ/rhoC+L/rhoC*df_s/dt
        dotQ = lambdaL*10 # NB better value?
        dfdT_Scheil = (f_s[i+1]-f_s[i])/dt
        print(dfdT_Scheil)
        T[i+1] = dT_next(dotQ,rhoC,L,dfdT_Scheil,dt,T[i])
        # --> Calc. new t_star due to temperature change
        t_st[i+1] = get_t_star(T_L[0],T[i+1],C_0[0],N,f_s[i+1],n_)
        
    t_plot = [i*dt for i in range (sim)]
    iteration = np.arange(sim)+1
    print('Scaled volume fraction X:')
    print(X)
    print('Volume fraction solid, alpha:')
    print(f_s)
    print('Temperature T:')
    print(T)
    print('t_star:')
    print(t_st)
    
    ##################################
    #plotting
    plt.figure(figsize=(6,6), dpi=120)
    plt.plot(t_plot,T,'r')
    plt.xlabel('t [s]')
    plt.ylabel(r'Temperature [$^{\circ}$C]')
    plt.title('Plot Temp')
    
    plt.figure(figsize=(6,6), dpi=120)
    plt.plot(t_plot,X,'b')
    plt.xlabel('t [s]')
    plt.ylabel(r'X')
    plt.title('Scaled volume fraction, X')

    plt.figure(figsize=(6,6), dpi=120)
    plt.plot(t_plot,f_s,'k')
    plt.xlabel('t [s]')
    plt.ylabel(r'f$_{s}$')
    plt.title(r'Volume fraction $\alpha$-phase, f$_{s}$')
    
    plt.figure(figsize=(6,6), dpi=120)
    plt.plot(iteration,t_st,'g')
    plt.xlabel('#')
    plt.ylabel(r't$_{star}$')
    plt.title('t$_{star}$ after each iteration')
    
    #subfig2[0].plot(t_plot,T,'b')
    #subfig2[1].plot(t_plot,X,'g')
    #subfig2[2].plot(t_plot,f_s,'k')
    #subfig2[3].plot(iteration,t_st,'r')
    
    # k = 1
    # for Xlist, n_id in zip(dXdt_anal,nlist):
    #     if k>3:
    #         subfig1[1].plot(t,Xlist, '--', label='n = {} (anal)'.format(n_id))
    #         k+=1
    #     else:
    #        subfig2[1].plot(t,Xlist, '--', label='n = {} (anal)'.format(n_id))
    #         k+=1
            
def main(argv):
    print('Program started')
    solidification()
    #subfig1[1].legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
    
#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
