#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 10:39:21 2022

@author: luiza.maciel
"""
from edos import *
from math import radians,degrees,sin,cos,sqrt,log,exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp,odeint
from constants import *

if __name__ == '__main__':
    
    V0=volume(teta0)
    T0=P0*10**-3*V0/(m_m*Rg)  
    Qa0=0
    Qp0=0
    W0=0
    
    x0=[P0,T0,Qa0,Qp0,W0] #condicoes iniciais
    print(x0)
    
    tk=np.linspace(teta0,tetaf,300)
    ti=tk[:find_nearest(tk,teta_comb)]
    tj=tk[find_nearest(tk,teta_comb):find_nearest(tk,teta_comb+delta_teta)]
    tf=tk[find_nearest(tk,teta_comb+delta_teta):]
    
    dt = (tetaf - teta0)/300
    
    P=np.ones(len(tk))*P0
    T=np.ones(len(tk))*T0
    Qa=np.ones(len(tk))*Qa0
    Qp=np.ones(len(tk))*Qp0
    W=np.ones(len(tk))*W0
    P1=np.ones(len(tk))*P0
    T1=np.ones(len(tk))*T
    
    compositions = np.zeros((len(tk),len(reactants_composition)))
    
################################################################################
#----------------------------------WIEBE------------------------------------#
    print('fcor= '+ str(fcor))
    for l in range(len(tk)-1):
        ts=[tk[l],tk[l+1]]
        x=solve_ivp(motored,ts,x0,args=(0,),method='DOP853').y
        
        P1[l+1]=x[0][-1]
        T1[l+1]=x[1][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    x0=[P0,T0,Qa0,Qp0,W0]
    #fechamento da valvula ao inicio da combustao
    for i in range(len(ti)-1): #range(0,98)
        ts=[ti[i],ti[i+1]]
        x=solve_ivp(motor,ts,x0,args=(0,),method='DOP853').y
        
        P[i+1]=x[0][-1] #ultima linha,coluna1
        T[i+1]=x[1][-1]
        Qa[i+1]=x[2][-1]
        Qp[i+1]=x[3][-1]
        W[i+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]

        
    for j in range (len(ti)-1,len(ti)+len(tj)-1):
        ts=[tk[j],tk[j+1]]
        
        xt = wiebe(tk[j])
        dxdt = dx(tk[j])
        
        x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt),method='DOP853').y

        P[j+1]=x[0][-1]
        T[j+1]=x[1][-1]
        Qa[j+1]=x[2][-1]
        Qp[j+1]=x[3][-1]
        W[j+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        massFraction.append(xt)
    
    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
        ts=[tk[f],tk[f+1]]
        
        xt = wiebe(tk[f])
        dxdt = dx(tk[f])
        
        x=solve_ivp(motorComb,ts,x0,args=(P1[f],xt,dxdt),method='DOP853').y
        
        P[f+1]=x[0][-1]
        T[f+1]=x[1][-1]
        Qa[f+1]=x[2][-1]
        Qp[f+1]=x[3][-1]
        W[f+1]=x[4][-1]
        
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        massFraction.append(xt)

        
        
################################################################################
#----------------------------------CINETICA------------------------------------#
#        
#    Spark = signal.unit_impulse(len(tk),find_nearest(tk,radians(-avanco)))*3000
#    
#    #curvaa de pressao sem combustao
#    for l in range(len(tk)-1):
#        ts=[tk[l],tk[l+1]]
#        spark=Spark[l]
#        x=solve_ivp(motored,ts,x0,args=(spark,),method='DOP853').y
#        
#        P1[l+1]=x[0][-1]
#        T1[l+1]=x[1][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#    
#    x0=[P0,T0,Qa0,Qp0,W0]
#    #fechamento da valvula ao inicio da combustao
#    for i in range(len(ti)-1): #range(0,98)
#        ts=[ti[i],ti[i+1]]
#        spark=Spark[i]
#        x=solve_ivp(motor,ts,x0,args=(spark,),method='DOP853').y
#        
#        P[i+1]=x[0][-1] #ultima linha,coluna1
#        T[i+1]=x[1][-1]
#        Qa[i+1]=x[2][-1]
#        Qp[i+1]=x[3][-1]
#        W[i+1]=x[4][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#        
#    #concentracao inicial - definidas aqui pois precisa de P e T
#    ro_fuel = P[i+1]* 10**-3 * M_c /(Ru* T[i+1])
#    ro_o2 = P[i+1]* 10**-3 * 32 /(Ru * T[i+1])
#    ro =  P[i+1]* 10**-3 * M_m /(Ru* T[i+1])
#    
#    c0fuel=(m_c/m_m)*ro*10**-3/M_c
#    c0O2=((m_air*0.232)/m_m)*ro*10**-3/32
#    concentration0 = [c0fuel,c0O2,0,0,0]
#    
#    xt=(Y0 - c0fuel * M_c *10**3/ro)/Y0
#    print(xt)
#    wiebeEv=[]
#    spark=0
#    #dxdt=(-dc(concentration0,ts[-1],T[i+1])[0]*M_c*10**3 / ro)/Y0
#    dxdt=0
##    #xt, dxdt = wiebe(ts[-1]), dx(ts[-1])
#    for j in range (len(ti)-1,len(ti)+len(tj)-1):
#        ts=[tk[j],tk[j+1]]
#        #spark = Spark[j]
#        x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt),method='DOP853').y
#        
#        P[j+1]=x[0][-1]
#        T[j+1]=x[1][-1]
#        Qa[j+1]=x[2][-1]
#        Qp[j+1]=x[3][-1]
#        W[j+1]=x[4][-1]
#        
#        xt,dxdt,next_c=twoStepKinects(concentration0,ts,T[j+1],P[j+1])
#        concentration0=next_c
#        print(concentration0)
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#       
#        massFraction.append(xt)
#    
#    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
#        ts=[tk[f],tk[f+1]]
#        #spark=Spark[f]
#        x=solve_ivp(motorComb,ts,x0,args=(P1[f],xt,dxdt),method='DOP853').y
#        
#        P[f+1]=x[0][-1]
#        T[f+1]=x[1][-1]
#        Qa[f+1]=x[2][-1]
#        Qp[f+1]=x[3][-1]
#        W[f+1]=x[4][-1]
#        
#        xt,dxdt,next_c=twoStepKinects(concentration0,ts,T[f+1],P[f+1])
#        concentration0=next_c
#
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#        
#        massFraction.append(xt)
#        
        
        
        
        
        
        
        
################################################################################
#----------------------------------EQUILIBRIO------------------------------------#
        
    #curvaa de pressao sem combustao
#    for l in range(len(tk)-1):
#        ts=[tk[l],tk[l+1]]
#        x=solve_ivp(motored,ts,x0,args=(0,),method='DOP853').y
#        
#        P1[l+1]=x[0][-1]
#        T1[l+1]=x[1][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#    
#    x0=[P0,T0,Qa0,Qp0,W0]
#    #fechamento da valvula ao inicio da combustao
#    for i in range(len(ti)-1): #range(0,98)
#        ts=[ti[i],ti[i+1]]
#        x=solve_ivp(motor,ts,x0,args=(0,),method='DOP853').y
#        
#        P[i+1]=x[0][-1] #ultima linha,coluna1
#        T[i+1]=x[1][-1]
#        Qa[i+1]=x[2][-1]
#        Qp[i+1]=x[3][-1]
#        W[i+1]=x[4][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#
#    next_composition=equilibrium(estimative,T[i+1],P[i+1]*9.869*10**-6)
#    composition0 = reactants_composition
#    massFraction.append(composition0[0])
#    
#    for j in range (len(ti)-1,len(ti)+len(tj)-1):
#        ts=[tk[j],tk[j+1]]
#    
#        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[j],composition0,next_composition,dt),method='DOP853').y
#        print(composition0,next_composition)
#        
#        P[j+1]=x[0][-1]
#        T[j+1]=x[1][-1]
#        Qa[j+1]=x[2][-1]
#        Qp[j+1]=x[3][-1]
#        W[j+1]=x[4][-1]
#        
#        composition0 = next_composition
#        estimative = composition0
#        
#        next_composition = equilibrium(estimative,T[j+1],P[j+1]*9.869*10**-6)
#        #print(next_composition)
#        massFraction.append(composition0[0])
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#        
#    
#    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
#        ts=[tk[f],tk[f+1]]
#        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[f],composition0,next_composition,dt),method='DOP853').y
#        
#        P[f+1]=x[0][-1]
#        T[f+1]=x[1][-1]
#        Qa[f+1]=x[2][-1]
#        Qp[f+1]=x[3][-1]
#        W[f+1]=x[4][-1]
#        
#        composition0 = next_composition
#        estimative = composition0
#        
#        next_composition = equilibrium(estimative,T[f+1],P[f+1]*9.869*10**-6)
#        massFraction.append(composition0[0])
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
    
################################################################################
#----------------------------------EQUILIBRIO 2D------------------------------------#
        
    #curvaa de pressao sem combustao
#    for l in range(len(tk)-1):
#        ts=[tk[l],tk[l+1]]
#        x=solve_ivp(motored,ts,x0,args=(0,),method='DOP853').y
#        
#        P1[l+1]=x[0][-1]
#        T1[l+1]=x[1][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#    
#    x0=[P0,T0,Qa0,Qp0,W0]
#    composition0 = reactants_composition
#    #fechamento da valvula ao inicio da combustao
#    for i in range(len(ti)-1): #range(0,98)
#        ts=[ti[i],ti[i+1]]
#        x=solve_ivp(motor,ts,x0,args=(0,),method='DOP853').y
#        
#        P[i+1]=x[0][-1] #ultima linha,coluna1
#        T[i+1]=x[1][-1]
#        Qa[i+1]=x[2][-1]
#        Qp[i+1]=x[3][-1]
#        W[i+1]=x[4][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#        
#        compositions[i] = reactants_composition
#   
#   
#    for j in range (len(ti)-1,len(ti)+len(tj)-1):
#        ts=[tk[j],tk[j+1]]
#    
#        burnedComposition = equilibrioAdiabatico(estimative, T[j],P[j]*9.869*10**-6,composition0,T[j])
#        
#        next_composition = wiebe(tk[j])*burnedComposition[1] + (1-wiebe(tk[j]))*reactants_composition
#        
#        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[j],composition0,next_composition,dt),method='DOP853').y
#        
#        
#        P[j+1]=x[0][-1]
#        T[j+1]=x[1][-1]
#        Qa[j+1]=x[2][-1]
#        Qp[j+1]=x[3][-1]
#        W[j+1]=x[4][-1]
#        
#        composition0 = next_composition
#        estimative = composition0
#        
#        massFraction.append(composition0[0])
#        compositions[j] = burnedComposition[2]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#        
#    
#    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
#        ts=[tk[f],tk[f+1]]
#        
#        burnedComposition = equilibrioAdiabatico(estimative, T[f],P[f]*9.869*10**-6,composition0, T[f])
#        next_composition = wiebe(tk[f])*burnedComposition[1] + (1-wiebe(tk[f]))*reactants_composition
#        
#        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[f],composition0,next_composition,dt),method='DOP853').y
#        
#        
#        P[f+1]=x[0][-1]
#        T[f+1]=x[1][-1]
#        Qa[f+1]=x[2][-1]
#        Qp[f+1]=x[3][-1]
#        W[f+1]=x[4][-1]
#        
#        composition0 = next_composition
#        estimative = composition0
#        
#        massFraction.append(composition0[0])
#        compositions[f] = burnedComposition[2]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]




################################################################################
#----------------------------GRAFICOS EQUILIBRIO-----------------------------#
    
    
    composicao_comb=compositions[:,0]
    composicao_o2=compositions[:,1]
    composicao_h2o=compositions[:,2]
    composicao_co2=compositions[:,3]
    composicao_co=compositions[:,4]
    composicao_o=compositions[:,5]
    composicao_n=compositions[:,6]
    composicao_no2=compositions[:,7]
    composicao_no=compositions[:,8]
    composicao_oh=compositions[:,9] 
    composicao_h2=compositions[:,10]
    composicao_h=compositions[:,11]
    composicao_n2=compositions[:,12]
    
    tkd=list(map(lambda x:degrees(x),tk))
    
    fig2=plt.figure(figsize=(15,10))
    
    plot=fig2.add_subplot(3,4,1)
    plot.plot(tkd,composicao_comb,'o')
    plot.legend(['Combustível'])
    
    plot=fig2.add_subplot(3,4,2)
    plot.plot(tkd,composicao_o2,'o')
    plot.legend(['O2'])
    
    plot=fig2.add_subplot(3,4,3)
    plot.plot(tkd,composicao_h2o,'o')
    plot.legend(['H2O'])
    
    plot=fig2.add_subplot(3,4,4)
    plot.plot(tkd,composicao_co2,'o')
    plot.legend(['CO2'])
    
    plot=fig2.add_subplot(3,4,5)
    plot.plot(tkd,composicao_co,'o')
    plot.legend(['CO'])
    
    plot=fig2.add_subplot(3,4,6)
    plot.plot(tkd,composicao_o,'o')
    plot.legend(['O'])
    
    plot=fig2.add_subplot(3,4,7)
    plot.plot(tkd,composicao_n,'o')
    plot.legend(['N'])
    
    plot=fig2.add_subplot(3,4,8)
    plot.plot(tkd,composicao_no2,'o')
    plot.legend(['NO2'])
    
    plot=fig2.add_subplot(3,4,9)
    plot.plot(tkd,composicao_no,'o')
    plot.legend(['NO'])
    
    plot=fig2.add_subplot(3,4,10)
    plot.plot(tkd,composicao_oh,'o')
    plot.legend(['OH'])
    
    plot=fig2.add_subplot(3,4,11)
    plot.plot(tkd,composicao_h2,'o')
    plot.legend(['H2'])
    
    plot=fig2.add_subplot(3,4,12)
    plot.plot(tkd,composicao_h,'o')
    plot.legend(['H'])
    
    plt.show()
    
    

###plots de graficos
    #print(P,T,Qa,Qp,W)
    print('\n'+'Pmax: '+ str(max(P)*10**-6)+' MPa' )
    print('Tmax: '+str(max(T))+' K' +str(degrees(tk[list(T).index(max(T))])))
    
    fig,(ax1,ax2)=plt.subplots(2,3)
    ax1[0].plot(tkd,P*10**(-6))
    ax1[1].plot(tkd,T)
    ax1[2].plot(tkd[len(ti)-1:-1],massFraction)
    ax2[0].plot(tkd,Qa)
    ax2[1].plot(tkd,Qp)
    ax2[2].plot(tkd,W)
    ax1[0].set_ylabel('P[MPa]')
    ax1[1].set_ylabel('T[K]')
    ax1[2].set_ylabel('X(teta)')
    ax2[0].set_ylabel('Qa[J]')
    ax2[1].set_ylabel('Qp[J]')
    ax2[2].set_ylabel('W[J]')
    plt.show()
