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
    massFraction=[0]
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
   
#    c0fuel=P[i]*X0/(Ru*T[i])
#    c0O2=P[i]*X0*OC/(Ru*T[i])
#    c0N2=P[i]*X0*OC*3.76/(Ru*T[i])
#    concentration0 = [c0fuel*10**-6,c0O2*10**-6,0,0,0,c0N2*10**-6] #kmol/cm3
   
    c0fuel=mols_c/volume(tk[i])
    c0O2=mols_c/volume(tk[i])*OC
    c0N2=mols_c/volume(tk[i])*OC*3.76
    
    concentration0 = [c0fuel*10**-3,c0O2*10**-3,0,0,0,c0N2*10**-3]
    
    for j in range (len(ti)-1,len(tk)-1):
        ts=[tk[j],tk[j+1]]
        spark = Spark[j]

       # xt,dxdt,next_c,molar_comp = twoStepKinects(concentration0,ts,T[j],P[j])
        xt,dxdt,next_c = twoStepKinectsVolume(concentration0,ts,T[j],P[j])
        print(xt,dxdt)
        
        x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt,spark),method='DOP853').y
      
        P[j+1]=x[0][-1]
        T[j+1]=x[1][-1]
        Qa[j+1]=x[2][-1]
        Qp[j+1]=x[3][-1]
        W[j+1]=x[4][-1]
        
        concentration0=next_c

        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
       
        massFraction.append(xt)
              
    return P,T,Qa,Qp,W,massFraction  