#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 15:05:16 2022

@author: luiza.maciel
"""
from edos import *
from math import radians,degrees,sin,cos,sqrt,log,exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp,odeint
from constants import *

################################################################################
#----------------------------------CINETICA------------------------------------#
def ODEkinects(ti,tj,tf,tk,x0,spark,tetaS,dt):
    
    P=np.ones(len(tk))*x0[0]
    T=np.ones(len(tk))*x0[1]
    Qa=np.ones(len(tk))*x0[2]
    Qp=np.ones(len(tk))*x0[3]
    W=np.ones(len(tk))*x0[4]
    
    P1=np.ones(len(tk))*x0[0]
    T1=np.ones(len(tk))*x0[1]
    
    P0,T0,Qa0,Qp0,W0=x0
    
    Spark = signal.unit_impulse(len(tk),find_nearest(tk,tetaS))*spark
    
    #curvaa de pressao sem combustao
    for l in range(len(tk)-1):
        ts=[tk[l],tk[l+1]]
        spark=Spark[l]
        x=solve_ivp(motored,ts,x0,method='DOP853').y
        
        P1[l+1]=x[0][-1]
        T1[l+1]=x[1][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
    # fechamento da valvula ao inicio da combustao  
    x0 = [P0,T0,Qa0,Qp0,W0]
    xt=0
    for i in range(len(ti)-1): #range(0,98)
        ts=[tk[i],tk[i+1]]
        spark=Spark[i]
        x=solve_ivp(motor,ts,x0,args=(spark,),method='DOP853').y
        
        P[i+1]=x[0][-1] #ultima linha,coluna1
        T[i+1]=x[1][-1]
        Qa[i+1]=x[2][-1]
        Qp[i+1]=x[3][-1]
        W[i+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        massFraction.append(xt)


   #concentracao inicial mol/m3 
    c0fuel=P[i]*X0/(Ru*T[i])
    c0O2=P[i]*X0*OC/(Ru*T[i])
    c0N2=P[i]*X0*OC*3.76/(Ru*T[i])
    concentration0 = [c0fuel*10**-6,c0O2*10**-6,0,0,0,c0N2*10**-6] #kmol/cm3
    cmp0=reactants_composition
    xt=(Y0 -  (c0fuel*Ru*T[i]/(P[i]))*(M_c/M_m))/Y0
    
    xt0=xt
    print(xt)
    dxdt=0

    for j in range (len(ti)-1,len(ti)+len(tj)-1):
        ts=[tk[j],tk[j+1]]
        spark = Spark[j]
        
        xt,dxdt,next_c,molar_comp = twoStepKinects(concentration0,ts,T[j],P[j])
        print(xt,dxdt)
        print(molar_comp)
        if xt>0.5 and xt0>xt or xt>1.0:
            xt,dxdt = xt0,0
            break

        cmp1 = [molar_comp[0],molar_comp[1],molar_comp[4],molar_comp[3],molar_comp[2],0,0,0,0,0,0,0,molar_comp[5]]
        #x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[j],cmp0,cmp1, dt),method='DOP853').y
        x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt,spark),method='DOP853').y
        xt0,dxdt0=xt,dxdt
        P[j+1]=x[0][-1]
        T[j+1]=x[1][-1]
        Qa[j+1]=x[2][-1]
        Qp[j+1]=x[3][-1]
        W[j+1]=x[4][-1]
        
        concentration0=next_c
        cmp0=cmp1
        xt0,dxdt0=xt,dxdt

        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
       
        massFraction.append(xt)
        
    print(degrees(tk[j]))
    for f in range(j,len(tk)-1): #200 ao 298
        ts=[tk[f],tk[f+1]]
        spark=Spark[f]
        
        x=solve_ivp(motorComb,ts,x0,args=(P1[f],xt,dxdt,spark),method='DOP853').y
        
        P[f+1]=x[0][-1]
        T[f+1]=x[1][-1]
        Qa[f+1]=x[2][-1]
        Qp[f+1]=x[3][-1]
        W[f+1]=x[4][-1]
        
        #xt,dxdt,next_c=twoStepKinects(concentration0,ts,T[f+1],P[f+1])
        #concentration0=next_c
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        massFraction.append(xt)
        
    return P,T,Qa,Qp,W,massFraction  
#    Spark = signal.unit_impulse(len(tk),find_nearest(tk,radians(-avanco)))*0
#    
#    #curvaa de pressao sem combustao
#    for l in range(len(tk)-1):
#        ts=[tk[l],tk[l+1]]
#        spark=Spark[l]
#        x=solve_ivp(motored,ts,x0,method='DOP853').y
#        
#        P1[l+1]=x[0][-1]
#        T1[l+1]=x[1][-1]
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#    
#    x0=[P0,T0,Qa0,Qp0,W0]
#   # fechamento da valvula ao inicio da combustao
#    xt=0
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
#        massFraction.append(xt)
#        
#    #concentracao inicial - definidas aqui pois precisa de P e T
#    
#    c0fuel=P[i]*X0/(Ru*T[i])
#    c0O2=P[i]*X0*OC/(Ru*T[i])
#    c0N2=P[i]*X0*OC*3.76/(Ru*T[i])
#    concentration0 = [c0fuel*10**6,c0O2*10**6,0,0,0]
#    composition0 = reactants_composition
#    
#    #xt=(Y0 -  (c0fuel*Ru*T0/P0)*(M_c/M_m))/Y0
#    xt=0
#    xt0=xt
#    print(xt)
#    wiebeEv=[]
#    spark=0
#    #dxdt=(-dc(concentration0,ts[-1],T[i+1])[0]*M_c*10**3 / ro)/Y0
#    dxdt=0
#
#    
#    for j in range (len(ti)-1,len(ti)+len(tj)-1):
#        ts=[tk[j],tk[j+1]]
#        spark = Spark[j]
#        
#        xt,dxdt,next_c = twoStepKinects(concentration0,ts,T[j],P[j])
#        if xt<xt0:
#            xt,dxdt=xt0,dxdt0
#        # dC,dO2,dCO,dCO2,dH2O
#        #h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2
#        c0N2=P[j]*X0*OC*3.76/(Ru*T[j])
#        cst0 = [concentration0[0],concentration0[1],concentration0[4],concentration0[3],concentration0[2],0,0,0,0,0,0,0,c0N2*10**6]
#        cst1 = [next_c[0],next_c[1],next_c[4],next_c[3],next_c[2],0,0,0,0,0,0,0,c0N2*10**6]
#        
#        
#        x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt),method='DOP853').y
#        
#        P[j+1]=x[0][-1]
#        T[j+1]=x[1][-1]
#        Qa[j+1]=x[2][-1]
#        Qp[j+1]=x[3][-1]
#        W[j+1]=x[4][-1]
#        
#        concentration0=next_c
#        xt0,dxdt0=xt,dxdt
#        print(xt,dxdt)
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#       
#        massFraction.append(xt)
#    
#    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
#        ts=[tk[f],tk[f+1]]
#        spark=Spark[f]
#        x=solve_ivp(motorComb,ts,x0,args=(P1[f],xt,dxdt),method='DOP853').y
#        
#        P[f+1]=x[0][-1]
#        T[f+1]=x[1][-1]
#        Qa[f+1]=x[2][-1]
#        Qp[f+1]=x[3][-1]
#        W[f+1]=x[4][-1]
#        
#        #xt,dxdt,next_c=twoStepKinects(concentration0,ts,T[f+1],P[f+1])
#        #concentration0=next_c
#
#        
#        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
#        
#        massFraction.append(xt)
#        return