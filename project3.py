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
    #listOfInput is orginaised as follows: [X_c, C_0, C_0_r, T_0, t_r, Nfrac, a, n, Scheil]
    #Where Nfrac is Nr/N
    #Both C_0 and N_frac might be list, but not at the same time. Then write [C_0_0, C_0_1, C_0_2...] in stead of one number. Will then only plot the temperature evolution

# Change of input parameters keeping all but one fixed:
listOfInput = [0.05, 1.1, 4.0, 670, 6, 1, 1, 3, False]#True]#False]                                    # <--- Standard input parameters
#listOfInput = [0.05, [2.0,3.2,4.8,6.0,7.2], 4.0, 670, 6, 1, 1, 3, True]                  # <--- Variation of the C_0/C_0_r ratio
#listOfInput = [0.05, 2.0, 4.0, 670, 6, [0.01,0.1,1,2,4,6,8,10,14], 1, 3, True]           # <--- Variation of the N_r/N ratio
#listOfInput = [0.05, 2.0, 4.0, 670, 6, 1, [0.8,1,1.2,1.3,1.4], 3, True]                  # <--- Variation of the external cooling rate a = L/(rho*c) <--- Low rates
#listOfInput = [0.05, 2.0, 4.0, 670, 6, 1, 1, [1,2,3], True]                              # <--- Variation of the time exponent n in the JMA-eq.

#PL = [dfdtlist0,Tlist1,Xlist2,flist3,fmlist4,Clist5,dTdtlist6] #PlotList
PI_glob = 1


#Global variables:
T_k = 273.15            #Degrees Kelvin at 0 degrees Centigrade
T_m = 660.5#+T_k         #Melting temperature of pure Al [K]
T_e = 577.0#+T_k         #Eutectic temperature of binary Al-Si [K]
C_s = 1.5               #Solubility of Si in Al at T_e, wt% Si
C_e = 12.2              #Eutectic composition, wt% Si


k_pc = C_s/C_e          #Partitioning coefficient defined to be C_sol/C_liq, is constant due to linearised phase diagram
m_upper = (T_m-T_e)/C_e #rate of linear line sol-liq
m_lower = (T_m-T_e)/C_s #rate of linear line sol-sol
#L = 1.0746				#Latent heat [J/mm^3]
L = 0.8                #Latent heat [J/mm^3]
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
t_sim = 250.0					#6 seconds simulation
dt = 0.01
Nt = math.ceil(t_sim/dt)
NiS = round(20/dt)


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
def getC_Scheil(C_0_t, f_m_t):
    return C_0_t*(1-f_m_t)**(k_pc-1)
def getC_Equi(C_0_t, f_m_t):
    return C_0_t/(1+(1-k_pc)*f_m_t)
#The reference temperature
def getT_r(T_L_t):
	return T_L_t - 2.0
#Time constant t_star (Applying the Hunt model, i.e. growth rate V prop. to undercooling^2/C_0
#NB: Define all parameters (calc from chosen reference condition)
def get_t_star(t_r_t, DT_r_t, DT_t, C_0_t, C_0_r_t, N_frac_t, f_m_t, f_m_r_t, n_t):
    return t_r_t*(DT_r_t/DT_t)**2*(C_0_t/C_0_r_t)*(N_frac_t)**(1/n_t)*(f_m_t/f_m_r_t)**(1/n_t)
    
def get_SF_eutectic(C_0_t):
    return C_e*math.exp(1/(k_pc-1))/C_0_t
 

######################## Solid weight fraction below #################
#The solid weight fraction at equilibrium, given by the lever rule.
def SF_Equi(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_t<T_S_t: return False #All is solid, do not distinguish between diferent solid phases
    return 1/(1-k_pc)*(T_L_t-T_t)/(T_m-T_t)
#The solid weight fraction differentiated with respect to temperature
def SF_Equi_dT(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_t<T_S_t: return 0 #All is solid, do not distinguish between diferent solid phases
    return 1/(k_pc-1)*(T_m-T_L_t)/(T_m-T_t)**2

#The solid weight fraction using the Scheil model. Non-equilibrium solidification.
def SF_Scheil(T_L_t, T_S_t, T_t):
    if T_L_t <= T_t:return 0
    if T_m < T_L_t:print('Error: T_L > T_melt'), exit()
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))> 1:print('Error: solid fraction > 1'), exit()
    if (1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1)))< 0:print('Error: negative solid fraction'), exit()
    return 1-((T_m-T_t)/(T_m-T_L_t))**(1/(k_pc-1))
#The solid weight fraction differentiated with respect to temperature
def SF_Scheil_dT(T_L_t, T_S_t, T_t):
    if T_t>T_L_t: return 0 #No solid has been formed at this temperature
    if T_m == T_L_t: return 1
    return (1/(k_pc-1))*(1/(T_m-T_L_t))*((T_m-T_t)/(T_m-T_L_t))**((2-k_pc)/(k_pc-1))

###################### Volume fraction ##########################
	
#Volume fraction as a function of a referance volume fraction and time, X_c and t_s respectively, time t and the integer n.
def XVF(X_c_t, n_t, t_t, t_s_t):
    return 1-(1-X_c_t)**((t_t/t_s_t)**n_t)
#Differentiate numerically XVF as a function of time, by 1st order finite difference method.
def dXVFdt(X_c_plus_t, X_c_minus_t, dt_t):
    return (X_c_plus_t-X_c_minus_t)/(2*dt_t)
#An analytic solution for dXVF/dt
def dXVFdt_anal(X_t,X_c_t,n_t,t_s_t):
    if X_t==0: return 0
    return -(1-X_t)*math.log(1-X_t)*n_t/(t_s_t*(math.log(1-X_t)/math.log(1-X_c_t))**(1/n_t))
#Precursor to differential form of MoleFraction
#Kick-off equation
def dXVFdt_precursor(n_t, t_s_t, X_c_t):
    return -(1-X_c_t)**((dt/t_s_t)**n_t)*math.log(1-X_c_t)*((dt/t_s_t)**n_t)*n_t/dt #Solve as backward Euler
############# Other functions ##################
#Total solidification time, eq. 15
def get_sol_time(t_eut_t,f_eut_t,L_t,a_t,rhoC_t):
    return t_eut_t+f_eut_t*L_t/(a_t*lambdaS/lambdaL)/rhoC_t
    
#Intrinsic cooling rate during solidification
def dT_next(dotQ_temp,rhoC_temp,L_temp,dfdT_Scheil_temp,dt_temp,T_prev_temp):
    return T_prev_temp-dt_temp*dotQ_temp/rhoC_temp+dt_temp*L_temp/rhoC_temp*dfdT_Scheil_temp

#Steady state cooling rate, i.e. dynamic balance occuring when X --> 1
def dT_next_steady_state(a_t,dfdT_Scheil_t):
    return a_t*lambdaS/lambdaL*(L/rhoC*dfdT_Scheil_t-1)**(-1)

#Weight function for smooth transition between transient and steady state evolution, centered around i_0. Width, alpha, is determined by numbers of points, NiS
def WF_smooth(i_t, i_0_t):
    alpha = NiS/16
    return(1/(1+math.exp(-(i_t-i_0_t)/alpha)))
    

def solidification(X_c, C_0, C_0_r, T_0, t_r, N_frac, a, n, Scheil, testPara = False, paraName='Default', PI = 1):
    if Scheil:
        SF_fun = SF_Scheil
        SF_fun_dT = SF_Scheil_dT
        getC_fun = getC_Scheil
    else:
        SF_fun = SF_Equi
        SF_fun_dT = SF_Equi_dT
        getC_fun = getC_Equi
############## Kickoff and initiation of itterations ##############

    #Must calculate variables which depend on reference parameters
    T_L_r = getT_L(C_0_r)
    T_S_r = getT_S(C_0_r) #This is used as a control in the solid weight fraction function (Scheil/Equilibrium)
    T_r = getT_r(T_L_r)
    f_m_r = SF_fun(T_L_r, T_S_r, T_r)
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
    f_m_0 = SF_fun(T_L, T_S, T_now)
    t_s_0 = get_t_star(t_r, DT_r, DT_0, C_0, C_0_r, N_frac, f_m_0, f_m_r, n)
    dXdt_now = dXVFdt_precursor(n, t_s_0, X_c)
    
    #Next temp (use variables for this time, i.e. t_0. No nucleation. i=0
    C_now = C_0
    f_m_now = SF_fun(T_L, T_S, T_now)
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

    #Initiate
    X_now = 0
    f_s_now = 0
    itt = round(t_sim/dt)
    bool_RSS = False #Reached steady state?

################ Starting iterations ###############
    for i in range(1,Nt):
        #Update from last iteration
        T_now = T_next      # T_now is for time j, while T_next is (j+1)
        
        #The following variables are defined from prev time iteration
        f_s_now = f_s_now+f_m_prev*dXdt_prev*dt
        if i == 1:      #since f_s_prev=0. Need to use kick off
           X_now = X_now + dXdt_prev*dt
        else:
            X_now = f_s_prev/f_m_prev
        C_now = getC_fun(C_0, f_m_prev)
        
        #These are defined from this time iteration
        DT_now  = T_L-T_now
        f_m_now = SF_fun(T_L, T_S, T_now)
        t_s_now = get_t_star(t_r, DT_r, DT_now, C_0, C_0_r, N_frac, f_m_now, f_m_r, n)
        dXdt_now = dXVFdt_anal(X_now, X_c, n, t_s_now)
        T_next = T_now-dt*a+dt*L/rhoC*f_m_now*dXdt_now

        #Redefine _prev for next iteration
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
        
        ########## Code control ###########
        if T_now < T_e or T_now > T_L:
            print('Temperature is crazy')
            print(T_L, T_e, T_now, 'T_L, T_e and T at i=',i)
            exit()
            break
        if f_s_now>f_m_now:
            print(i,'Transient, f_s>f_m')
            exit()
            break
        if X_now > 1-1e-1 and not bool_RSS:
            bool_RSS = True
            itt = i
            print(t_0+i*dt,'Done with Transient, at i={}'.format(i),T_now)
            print(f_m_now, f_s_now, f_s_now/f_m_now)
        if T_now < T_e:
            print(i)
            break
        if bool_RSS and i-itt > NiS: #NiS times into steady state with transient, we smooth
            break

    ################# Time for steady state evolution of our system ##############
    DoneSmooth = False
    dT_prev = dTdtlist[itt]*dt
    dfmdT_prev = dfdtlist[itt]*dt/dT_prev
    f_s_prev = flist[itt]
    f_m_prev = fmlist[itt]
    T_next = Tlist[itt]
    BoolAbort = False
    if bool_RSS and 1:
        while T_next > T_e and timelist[-1] < t_sim:
            if itt == i+2:
                DoneSmooth = True
            T_now = T_next

            #Defined from prev
            X_now = f_s_prev/f_m_prev
            C_now = getC_fun(C_0,f_m_prev)
            dfm = dfmdT_prev*dT_prev
            f_s_now = f_s_prev+dfm
                
            #Update from this iteration
            f_m_now = SF_fun(T_L, T_S, T_now)
            dfmdT = SF_fun_dT(T_L, T_S, T_now)
            dTdt = dT_next_steady_state(a, dfmdT)
            dT_now = dTdt*dt

            #Update variables
            T_prev = T_now
            f_m_prev = f_m_now
            if not f_m_now:
                print('Temperature is below solidus temperature in the Lever rule method. We are now in a one phase section of the phase diagram.')
                BoolAbort = True
                Tlist = Tlist + [Tlist[-1],Tlist[-1]-(t_sim-timelist[-1])*a*lambdaS/lambdaL]
                timelist = timelist + [timelist[-1],t_sim]
                dfdtlist = dfdtlist + [0,0]
                flist = flist + [1,1]
                Xlist = Xlist + [1,1]
                Clist = Clist + [Clist[-1],Clist[-1]]
                fmlist = fmlist + [1,1]
                dTdtlist = dTdtlist + [a*lambdaS/lambdaL,a*lambdaS/lambdaL]

                break
            f_s_prev = f_s_now
            dT_prev = dT_now
            dfmdT_prev = dfmdT
            T_next = T_now + dT_now

            if not DoneSmooth:
                WF = WF_smooth((itt-i+NiS),round(NiS/2))
                C_now = (1-WF)*Clist[itt]+WF*C_now
                f_m_now = (1-WF)*fmlist[itt]+WF*f_m_now
                f_s_now = (1-WF)*flist[itt]+WF*f_s_now
                dTdt = (1-WF)*dTdtlist[itt]+WF*dTdt
                dfm = (1-WF)*dfdtlist[itt]*dt+WF*dfm
                T_now = (1-WF)*Tlist[itt] + WF*T_next
            #    print(i,itt, dTdt,dTdtlist[itt], WF)

            

            if DoneSmooth:
                Tlist.append(T_now)
                timelist.append(t_0+itt*dt)
                dfdtlist.append(dfm/dt) #This should be fixed to the previous step
                flist.append(f_s_now)
                Xlist.append(X_now)
                Clist.append(C_now)
                fmlist.append(f_m_now)
                dTdtlist.append(dTdt)
            else:
                Tlist[itt] = T_now
                timelist[itt] = t_0+itt*dt
                dfdtlist[itt] = dfm/dt #This should be fixed to the previous step
                flist[itt] = f_s_now
                Xlist[itt] = X_now
                Clist[itt] = C_now
                fmlist[itt] = f_m_now
                dTdtlist[itt] = dTdt
            #Keep track of time
            itt = itt+1
           

            ############### Code control ################
            if f_s_now > f_m_now:
                print(itt, 'f_s > f_m')
                break
            if T_now > T_L:
                print('T_now > T_L',itt)
                exit()
    ################# Eutectic point reached ##############
    # Here Gibb's phase rule predicts that the binary Al-Si eutectic must proceed at constant temperature, i.e. dT/dt = 0.
    #timelist.append(t_sim)
    #Tlist.append(T_e)
    if bool_RSS and timelist[-1] < t_sim and not BoolAbort:
        T_now = T_e

        f_eut = flist[-1]
        t_eut = get_sol_time(0,f_eut,L,a,rhoC)
        t_f = get_sol_time(timelist[-1],f_eut,L,a,rhoC)
        print('Duration of eutectic solidification: {} s'.format(t_eut))
        tfLTts = True # t_f Less Than t_sim
        if t_f > t_sim:
            tfLTts = False
            t_f = t_sim
        
        Tlist = Tlist + [T_e, T_e]
        timelist = timelist + [timelist[-1], t_f]
        dfdtlist = dfdtlist + [(1-flist[-1])/t_eut,(1-flist[-1])/t_eut]
        flist = flist + [flist[-1],1]
        print(flist[-1],flist[-2])
        Xlist = Xlist + [Xlist[-1], Xlist[-1]]
        Clist = Clist + [Clist[-1], Clist[-1]]
        fmlist = fmlist + [1,1]
        dTdtlist = dTdtlist + [0,0]
        
        if tfLTts:
            Tlist = Tlist + [T_e, T_e-(t_sim-t_f)*a*lambdaS/lambdaL]
            timelist = timelist + [t_f+dt, t_sim]
            dfdtlist = dfdtlist + [0,0]
            flist = flist + [1,1]
            Xlist = Xlist + [Xlist[-1], Xlist[-1]]
            Clist = Clist + [Clist[-1], Clist[-1]]
            fmlist = fmlist + [1,1]
            dTdtlist = dTdtlist + [a*lambdaS/lambdaL,a*lambdaS/lambdaL]
            

    ################    Plotting     ######################
    
    SF = 0 #Samefig, executes subplot which does not share yscale. For the T-dfdt plot
 #   PB = [1,1,1,1,1,1,1] #PlotBool
    PB = [0,1,0,1,1,0,0] #PlotBool
    PL = [dfdtlist,Tlist,Xlist,flist,fmlist,Clist,dTdtlist] #PlotList
    PN = ['Evolution of solidificationrate','Temperature evolution','Scaled volume fraction evolution',\
    'Evolution of volume fraction formed', 'Evolution of maximum theoretical volume fraction', 'Evolution of Si wt% concentration in liquid', 'Evolution of temperature gradient']   #PlotNames
    PY = ['df/dt',r'Temperature [$^{\circ}$C]','X', r'f$_{s}$', 'f$_{m}$', r'C$_{L}$ [wt% Si]', 'dT/dt [$^{\circ}$C/s]']
    PX = 't [s]'
    print('\n\n\n\nLength of vectors:\n')
    print('Timelist:')
    print(np.size(timelist))
    print('dfdtlist:')
    print(np.size(dfdtlist))
    if testPara:
        plt.plot(timelist, PL[PI], label = paraName)
        plt.title(PN[PI],fontsize= 30,y=1.04)
        plt.xlabel(PX)
        plt.ylabel(PY[PI])
        params = {'font.size': 28, 'legend.fontsize': 24,'lines.linewidth': 3.0}
        plt.rcParams.update(params)
        return 0 
    if SF:
        firstname = True
        for PB_t,PN_t in zip(PB,PN):
            if firstname and PB_t:
                origin = PN_t
                exec("fig, "+PN_t+" = plt.subplots()")
                firstname = False
            elif PB_t:
                exec(PN_t+"="+origin+".twinx()")
    for PB_t, PL_t, PN_t, PY_t in zip(PB,PL,PN, PY):
        if PB_t:
            if SF:
                exec(PN_t+".plot(timelist,PL_t)")
                exec(PN_t+".set_ylabel(PY_t)")
                exec(PN_t+".set_xlabel(PX)")
            else:
                plt.figure()
                plt.plot(timelist,PL_t)
                plt.title(PN_t,fontsize= 16,y=1.04)
                plt.xlabel(PX)
                plt.ylabel(PY_t)
                plt.rcParams.update({'font.size': 16})
def main(argv):
    #listOfInput is orginaised as follows: [X_c, C_0, C_0_r, T_0, t_r, Nfrac, a, n, Scheil]
    #PL = [dfdtlist,Tlist,Xlist,flist,fmlist,Clist,dTdtlist] #PlotList
    loi = listOfInput
    if loi[-1]:
        func = "Scheil method"
    else:
        func = "Lever-rule method"
    print('Program started. Using the ', func)
    if type(loi[1]) == type([]):
        plt.figure()
        for C_0_ in loi[1]:
            solidification(loi[0], C_0_, loi[2], loi[3], loi[4], loi[5], loi[6], loi[7], loi[8], True, r"C$_{0}$="+str(C_0_),PI_glob)
    elif type(loi[5]) == type([]):
        plt.figure()
        for N_frac_ in loi[5]:
            print('Started running using standard input parameters and N_frac={}\n'.format(N_frac_))
            solidification(loi[0], loi[1], loi[2], loi[3], loi[4], N_frac_, loi[6], loi[7], loi[8], True, r"N$_{r}$/N="+str(N_frac_),PI_glob)
    elif type(loi[6]) == type([]):
        plt.figure()
        for a in loi[6]:
            print('Started running using standard input parameters and a={}\n'.format(a))
            solidification(loi[0], loi[1], loi[2], loi[3], loi[4], loi[5], a, loi[7], loi[8], True, "a="+str(a),PI_glob)
    elif type(loi[7]) == type([]):
        plt.figure()
        for n in loi[7]:
            print('Started running using standard input parameters and n={}\n'.format(n))
            solidification(loi[0], loi[1], loi[2], loi[3], loi[4], loi[5], loi[6], n, loi[8], True, "n="+str(n),PI_glob)
    else:
        solidification(loi[0], loi[1], loi[2], loi[3], loi[4], loi[5], loi[6], loi[7], loi[8])
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    #plt.legend()
    plt.show()
    
#Only run if this is a main file, and not a module
if __name__ == "__main__":
    main(sys.argv[1:])
