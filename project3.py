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

#List of input parameters, will be initiated in the main function.
#Order of ...
listOfInput = [0.05, 2.0, 4.0]


#Global variables:
T_k = 273.15            #Degrees Kelvin at 0 degrees Centigrade
T_m = 660.5#+T_k         #Melting temperature of pure Al [K]
T_e = 577.0#+T_k         #Eutectic temperature of binary Al-Si [K]
C_s = 1.5               #Solubility of Si in Al at T_e, wt% Si
C_e = 12.2              #Eutectic composition, wt% Si


k_pc = C_s/C_e          #Partitioning coefficient defined to be C_sol/C_liq, is constant due to linearised phase diagram
m_upper = (T_m-T_e)/C_e #rate of linear line sol-liq
m_lower = (T_m-T_e)/C_s #rate of linear line sol-sol
L = 1.0746				#Latent heat [J/mm^3]
rho = 2.7e-3           #Density[g/mm^3]
C_hc = 24.20            #Heat capacity [J/mol]
M_mAl = 26.98          #[g/mol] Molar mass Aluminium
rhoC = rho/M_mAl*C_hc #Volume heat capacity [J/(Celcius*mm^3)]
#print('rhoC1 = {0:.5f} J/Cmm^3'.format(rhoC2))
rhoC2 = 0.0027			   #Volume heat capacity [J/(Celcius*mm^3)] density dependent
print('rhoC = {0:.5f} J/Cmm^3'.format(rhoC))
#Andre verdiar frå literaturen:
lambdaL = 0.094			#Thermal conductivity liquid Al [W/(Celcius*mm)] @ 665 degrees Celcius
lambdaS = 0.213			#Thermal conductivity solid Al [W/(Celcius*mm)] @ 630 degrees Celcius
t_sim = 320.0					#6 seconds simulation
dt = 0.01
Nt = math.ceil(t_sim/dt)

##################### All getStuff is defined here ##############
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
def getC_scheil(C_0_t, f_m_t):
    return C_0_t*(1-f_m_t)**(k_pc-1)
#The reference temperature
def getT_r(T_L_t):
	return T_L_t - 2.0
#Time constant t_star (Applying the Hunt model, i.e. growth rate V prop. to undercooling^2/C_0
#NB: Define all parameters (calc from chosen reference condition)
def get_t_star(t_r_t, DT_r_t, DT_t, C_0_t, C_0_r_t, N_r_t, N_t, f_m_t, f_m_r_t, n_t):
    return t_r_t*(DT_r_t/DT_t)**2*(C_0_t/C_0_r_t)*(N_r_t/N_t)**(1/n_t)*(f_m_t/f_m_r_t)**(1/n_t)
 

######################## Solid weight fraction below #################
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
    if T_L_t <= T_t:return 0
    if T_m < T_L_t:print('Error: T_L > T_melt'), exit()
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))> 1:print('Error: solid fraction > 1'), exit()
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))< 0:print('Error: negative solid fraction'), exit()
    return 1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1))
#The solid weight fraction differentiated with respect to temperature
def SF_scheil_dT(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_m == T_L_t: return 1
    return (1/(k_pc-1))*(1/(T_m-T_L_t))*((T_m-T_t)/(T_m-T_L_t))**((2-k_pc)/(k_pc-1))
def SF_eutectic(C_0_t):
    return C_e*math.exp(1/(k_pc-1))/C_0_t

###################### Volume fraction ##########################
	
#Volume fraction as a function of a referance volume fraction and time, X_c and t_s respectively, time t and the integer n.
def XVF(X_c_t, n, t_t, t_s_t):
    return 1-(1-X_c_t)**((t_t/t_s_t)**n)
#Differentiate numerically XVF as a function of time, by 1st order finite difference method.
def dXVFdt(X_c_plus_t, X_c_minus_t, dt_t):
    return (X_c_plus_t-X_c_minus_t)/(2*dt_t)
#An analytic solution for dXVF/dt
def dXVFdt_anal(X_t,X_c_t,n,t_s_t):
    if X_t==0: return 0
    return -(1-X_t)*math.log(1-X_t)*n/(t_s_t*(math.log(1-X_t)/math.log(1-X_c_t))**(1/n))
#Precursor to differential form of MoleFraction
#Kick-off equation
def dXVFdt_precursor(n_t, t_s_t, X_c_t):
    return -(1-X_c_t)**((dt/t_s_t)**n_t)*math.log(1-X_c_t)*((dt/t_s_t)**n_t)*n_t/dt #Solve as backward Euler
############# Other functions ##################
#Total solidification time, eq. 15
def sol_time(t_eq_temp,f_eq_temp,L_temp,a_max_temp,rhoC_temp):
    return t_eq_temp-f_eq_temp*L_temp/a_max_temp/rhoC_temp
    
#Intrinsic cooling rate during solidification
def dT_next(dotQ_temp,rhoC_temp,L_temp,dfdT_Scheil_temp,dt_temp,T_prev_temp):
    return T_prev_temp-dt_temp*dotQ_temp/rhoC_temp+dt_temp*L_temp/rhoC_temp*dfdT_Scheil_temp

#Steady state cooling rate, i.e. dynamic balance occuring when X --> 1
def dT_next_steady_state(a_t,dfdT_Scheil_t):
    return a_t*lambdaS/lambdaL*(L/rhoC*dfdT_Scheil_t-1)**(-1)

def solidification(X_c, C_0, C_0_r, T_0 = 670, t_r=6, N=1000, N_r=1000, a=1.0, n=3):
    #Must calculate variables which depend on reference parameters
    T_L_r = getT_L(C_0_r)
    T_S_r = getT_S(C_0_r) #This is used as a control in the solid weight fraction function (Scheil/Equilibrium)
    T_r = getT_r(T_L_r)
    f_m_r = SF_scheil(T_L_r, T_S_r, T_r)
    DT_r = T_L_r-T_r #Undercooling, is equal to 2 degrees
    
    #Constants of our system. From chosen initial values
    T_L = getT_L(C_0)
    T_S = getT_S(C_0)
    T_n = getT_n(T_L)
    t_n = (T_L-T_n)/a

    #Undercooling kickoff
    DT_0 = T_L-T_n #Is 0.1

    #kickoff i = 0
    T_now = T_n
    t_0 = (T_0-T_n)/a
    f_m_0 = SF_scheil(T_L, T_S, T_now)
 #   print((f_m_0/(2*f_m_r))**(1/3)), exit()
 #   print(f_m_now, 'f_m')
    t_s_0 = get_t_star(t_r, DT_r, DT_0, C_0, C_0_r, N_r, N, f_m_0, f_m_r, n)
    X_now = 0
    dXdt_ko = dXVFdt_precursor(n, t_s_0, X_c)
    dXdt_now = dXdt_ko
    
    #Next temp (use variables for this time, i.e. t_0. No nucleation. i=0
    C_now = C_0
    T_L_now = getT_L(C_now)
    T_S_now = getT_S(C_now)
    f_m_now = SF_scheil(T_L_now, T_S_now, T_now)
    T_next = T_now-dt*a+dt*L/rhoC*f_m_now*dXdt_now #i=1, from i=0 values

    #All X related should be determined from last step, we need prev values
    f_m_prev = f_m_now
    f_s_prev = 0
    dXdt_prev = dXdt_now
    T_prev = T_now

    #Different lists
    Tlist = [T_0,T_n]
    timelist =[0,t_0]
    dfdtlist = [0,0]
    flist = [0,0]
    Xlist = [0,0]
    Clist = [C_0, C_now]
    fmlist = [f_m_now, f_m_now]
    dTdtlist = [-a,-a]

    X_now = 0 #initiate
    f_s_now = 0
    itt = 0
    bool_RSS = False #Reached steady state?

################ Starting itterations ###############
    for i in range(1,Nt):
        
        #Update from last itteration
        T_now = T_next      # T_now is for time j, while T_next is (j+1)
        
        #The following variables are defined from prev time itteration
        f_s_now = f_s_now+f_m_prev*dXdt_prev*dt
        X_now = X_now + dXdt_prev*dt
        C_now = getC_scheil(C_0, f_m_prev)
        
        #These are defined from this time itteration
        DT_now  = T_L-T_now
        f_m_now = SF_scheil(T_L, T_S, T_now)
        t_s_now = get_t_star(t_r, DT_r, DT_now, C_0, C_0_r, N_r, N, f_m_now, f_m_r, n)
        dXdt_now = dXVFdt_anal(X_now, X_c, n, t_s_now)
        T_next = T_now-dt*a+dt*L/rhoC*f_m_now*dXdt_now

        #Redefine _prev for next itteration
        f_m_prev = f_m_now
        f_s_prev = f_s_now
        dXdt_prev = dXdt_now
        T_prev = T_now


        Tlist.append(T_now)
        timelist.append(t_0+i*dt)
        dfdtlist.append(f_m_now*dXdt_now)
        flist.append(f_s_now)
        Xlist.append(X_now)
        Clist.append(C_now)
        fmlist.append(f_m_now)
        dTdtlist.append((-a+L/rhoC*f_m_now*dXdt_now))
        
        #Code control
        if T_now < T_e or T_now > T_L:
            print('Temperature is crazy')
            print(T_L, T_e, T_now, 'T_L, T_e and T at i=',i)
            exit()
            break
        if f_s_now>f_m_now:
            print(i,'Transient, f_s>f_m')
            exit()
            break
        if X_now > 1-5e-2:
            bool_RSS = True
            itt = i
            print(t_0+i*dt,'Done with Transient, at i={}'.format(i),T_now)
            print(f_m_now, f_s_now, f_s_now/f_m_now)
            break
        if T_now < T_e:
            print(i)
            break
    if bool_RSS and 1:
        while T_now > T_e:
            C_now = getC_scheil(C_0,f_m_now)#C_0*(1-f_m_now)**(k_pc-1)
            f_m_now = SF_scheil(T_L, T_S, T_now)
            dfdT = SF_scheil_dT(T_L, T_S, T_now)
            dTdt = dT_next_steady_state(a, dfdT)
            dT = dTdt*dt
            df = dfdT*dTdt*dt
            f_s_now = f_s_now+df
            T_now = T_now + dT
            X_now = f_s_now/f_m_now
            itt = itt+1
            Tlist.append(T_now)
            timelist.append(t_0+itt*dt)
            dfdtlist.append(df/dt) #This should be fixed to the previous step
            flist.append(f_s_now)
            Xlist.append(X_now)
            Clist.append(C_now)
            fmlist.append(f_m_now)
            dTdtlist.append(dTdt)
            if dfdtlist[-2]/dfdtlist[-1] < 1.5 and False:
                print('dfdt is not continous at all', itt)
                exit()
            if f_s_now > f_m_now:
                print(itt, 'f_s > f_m')
                break
            if T_now > T_L_now:
                print('dammit',itt)
                exit()
            if itt*dt > t_sim:
                print('Never reached T_e')
                break
    PB = [1,1,0,0,0,0,0] #PlotBool
    PB = [1,1,1,1,1,1,1] #PlotBool
    PL = [dfdtlist,Tlist,Xlist,flist,fmlist,Clist,dTdtlist] #PlotList
    PN = ['dfdt','Temp','X', 'f_s', 'f_m', 'C_L','dTdt'] #PlotNames
    SF = 0 #Samefig, executes subplot which does not share yscale. For the T-dfdt plot
    if SF:
        firstname = True
        for PB_t,PN_t in zip(PB,PN):
            if firstname and PB_t:
                origin = PN_t
                exec("fig, "+PN_t+" = plt.subplots()")
                firstname = False
            elif PB_t:
                exec(PN_t+"="+origin+".twinx()")
    for PB_t, PL_t, PN_t in zip(PB,PL,PN):
        if PB_t:
            if SF:
                exec(PN_t+".plot(timelist,PL_t)")
            else:
                plt.figure()
                plt.plot(timelist,PL_t)
                plt.title(PN_t)
                
    
    plt.show()

#################### Better solution above
    exit()

    
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
    
            
def main(argv):
    print('Program started')
    solidification(0.05, 2.0, 4.0)
    plt.show()
    
#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
