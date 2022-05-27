#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  7 20:37:40 2022

@author: luiza.maciel
"""

from edos import *
from math import radians,degrees,sin,cos,sqrt,log,exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp,odeint
from constants import *
from equilibrium import *

################################################################################
#----------------------------------EQUILIBRIO------------------------------------#
def ODEequilibrium1D(ti,tj,tf,tk,x0,reactants_composition,estimative,dt):     
    P=np.ones(len(tk))*x0[0]
    T=np.ones(len(tk))*x0[1]
    Qa=np.ones(len(tk))*x0[2]
    Qp=np.ones(len(tk))*x0[3]
    W=np.ones(len(tk))*x0[4]
    
    P1=np.ones(len(tk))*x0[0]
    T1=np.ones(len(tk))*x0[1]
    
    P0,T0,Qa0,Qp0,W0=x0
    compositions = np.zeros((len(tk),len(reactants_composition)))
    k_eval = []       
    
    #curvaa de pressao sem combustao
    for l in range(len(tk)-1):
        ts=[tk[l],tk[l+1]]
        x=solve_ivp(motored,ts,x0,method='DOP853').y
        
        P1[l+1]=x[0][-1]
        T1[l+1]=x[1][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    x0 = [P0,T0,Qa0,Qp0,W0]
    composition0 = reactants_composition
    #fechamento da valvula ao inicio da combustao
    for i in range(len(ti)-1): #range(0,98)
        ts=[ti[i],ti[i+1]]
        x=solve_ivp(motor,ts,x0,args=(0,),method='DOP853').y
        
        ke = K(T[i],composition0)[0]
        k_eval.append(ke)
        
        P[i+1]=x[0][-1] #ultima linha,coluna1
        T[i+1]=x[1][-1]
        Qa[i+1]=x[2][-1]
        Qp[i+1]=x[3][-1]
        W[i+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        compositions[i] = np.asarray(reactants_composition)/sum(reactants_composition)
   
   #combustao
    for j in range (len(ti)-1,len(ti)+len(tj)-1):
        ts=[tk[j],tk[j+1]]
        
        if dvdt(tk[j]) < 1e-5:
            burned_composition = equilibrioAdiabaticoVolumeConstante(estimative, T[j],P[j]*9.869*10**-6,volume(tk[j]),composition0,T[j])[1]
        else:
            burned_composition = equilibrioAdiabatico(estimative, T[j],P[j]*9.869*10**-6,composition0,T[j])[1]
       # burnedComposition = equilibrium(estimative, T[j],P[j]*9.869*10**-6)
        
        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[j],composition0,burned_composition,dt),method='DOP853').y
        
        ke = K(T[j],composition0)[0]
        k_eval.append(ke)
        
        P[j+1]=x[0][-1]
        T[j+1]=x[1][-1]
        Qa[j+1]=x[2][-1]
        Qp[j+1]=x[3][-1]
        W[j+1]=x[4][-1]
        print(burned_composition[0])
        composition0 = burned_composition
        estimative = composition0
        
        compositions[j] = np.asarray(burned_composition)/sum(burned_composition)
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
    
    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
        ts=[tk[f],tk[f+1]]
        
        if dvdt(tk[f]) < 1e-5:
            burned_composition = equilibrioAdiabaticoVolumeConstante(estimative, T[j],P[j]*9.869*10**-6,volume(tk[f]),composition0,T[f])[1]
        else:
            burned_composition = equilibrioAdiabatico(estimative, T[f],P[f]*9.869*10**-6,composition0,T[f])[1]
        
        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[f],composition0,burned_composition,dt),method='DOP853').y
        
        ke = K(T[f],composition0)[0]
        k_eval.append(ke)
        
        P[f+1]=x[0][-1]
        T[f+1]=x[1][-1]
        Qa[f+1]=x[2][-1]
        Qp[f+1]=x[3][-1]
        W[f+1]=x[4][-1]
        
        composition0 = burned_composition
        estimative = composition0
        
        compositions[f] = np.asarray(burned_composition)/sum(burned_composition)
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    return P,T,Qa,Qp,W,k_eval,compositions
        
    
###############################################################################
#----------------------------------EQUILIBRIO 2D------------------------------------#
def ODEequilibrium2D(ti,tj,tf,tk,x0,reactants_composition,estimative,dt):    
    P=np.ones(len(tk))*x0[0]
    T=np.ones(len(tk))*x0[1]
    Qa=np.ones(len(tk))*x0[2]
    Qp=np.ones(len(tk))*x0[3]
    W=np.ones(len(tk))*x0[4]
    
    P1=np.ones(len(tk))*x0[0]
    T1=np.ones(len(tk))*x0[1]
    
    P0,T0,Qa0,Qp0,W0=x0
    compositions = np.zeros((len(tk),len(reactants_composition)))
    k_eval = []       
    
    #curvaa de pressao sem combustao
    for l in range(len(tk)-1):
        ts=[tk[l],tk[l+1]]
        x=solve_ivp(motored,ts,x0,method='DOP853').y
        
        P1[l+1]=x[0][-1]
        T1[l+1]=x[1][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    x0 = [P0,T0,Qa0,Qp0,W0]
    composition0 = reactants_composition
    
    #fechamento da valvula ao inicio da combustao
    for i in range(len(ti)-1): #range(0,98)
        ts=[ti[i],ti[i+1]]
        x=solve_ivp(motor,ts,x0,args=(0,),method='DOP853').y
        
        ke = K(T[i],composition0)[0]
        k_eval.append(ke)
        
        P[i+1]=x[0][-1] #ultima linha,coluna1
        T[i+1]=x[1][-1]
        Qa[i+1]=x[2][-1]
        Qp[i+1]=x[3][-1]
        W[i+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        compositions[i] = np.asarray(reactants_composition)/sum(reactants_composition)
   
   #combustao
    for j in range (len(ti)-1,len(ti)+len(tj)-1):
        ts=[tk[j],tk[j+1]]
        
        
        if dvdt(tk[j]) < 1e-5:
            burnedComposition = equilibrioAdiabaticoVolumeConstante(estimative, T[j],P[j]*9.869*10**-6,volume(tk[j]),composition0,T[j])
        else:
            burnedComposition = equilibrioAdiabatico(estimative, T[j],P[j]*9.869*10**-6,composition0,T[j])
       # burnedComposition = equilibrium(estimative, T[j],P[j]*9.869*10**-6)
        next_composition = wiebe(tk[j])*burnedComposition[1] + (1-wiebe(tk[j]))*reactants_composition
        
        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[j],composition0,next_composition,dt),method='DOP853').y
        
        ke = K(T[j],composition0)[0]
        k_eval.append(ke)
        
        P[j+1]=x[0][-1]
        T[j+1]=x[1][-1]
        Qa[j+1]=x[2][-1]
        Qp[j+1]=x[3][-1]
        W[j+1]=x[4][-1]
        print(next_composition[0])
        composition0 = next_composition
        estimative = composition0
        
        compositions[j] = np.asarray(next_composition)/sum(next_composition)
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
    
    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
        ts=[tk[f],tk[f+1]]
        if dvdt(tk[f]) < 1e-5:
            burnedComposition = equilibrioAdiabaticoVolumeConstante(estimative, T[f],P[f]*9.869*10**-6,volume(tk[f]),composition0,T[f])
        else:
            burnedComposition = equilibrioAdiabatico(estimative, T[f],P[f]*9.869*10**-6,composition0,T[f])
    
        next_composition = wiebe(tk[f])*burnedComposition[1] + (1-wiebe(tk[f]))*reactants_composition
        
        x=solve_ivp(motorCombEquilibrium,ts,x0,args=(P1[f],composition0,next_composition,dt),method='DOP853').y
        
        ke = K(T[f],composition0)[0]
        k_eval.append(ke)
        
        P[f+1]=x[0][-1]
        T[f+1]=x[1][-1]
        Qa[f+1]=x[2][-1]
        Qp[f+1]=x[3][-1]
        W[f+1]=x[4][-1]
        
        composition0 = next_composition
        estimative = composition0
        
        compositions[f] = np.asarray(next_composition)/sum(next_composition)
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    
    return P,T,Qa,Qp,W,k_eval,compositions
