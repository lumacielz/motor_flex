#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  4 14:23:47 2022

@author: luiza.maciel
"""
from edos import *
from math import radians,degrees,sin,cos,sqrt,log,exp
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp,odeint
from constants import *

################################################################################
#----------------------------------WIEBE------------------------------------#
def ODEwiebe(ti,tj,tf,tk,x0):
   
    P=np.ones(len(tk))*x0[0]
    T=np.ones(len(tk))*x0[1]
    Qa=np.ones(len(tk))*x0[2]
    Qp=np.ones(len(tk))*x0[3]
    W=np.ones(len(tk))*x0[4]
    
    P1=np.ones(len(tk))*x0[0]
    T1=np.ones(len(tk))*x0[1]
    P0,T0,Qa0,Qp0,W0=x0
    
    massFraction=np.zeros(len(tk))
    k_eval = []
    xt=0

    #calculo da curva de pressao sem combustao
    for l in range(len(tk)-1):
        ts=[tk[l],tk[l+1]]
        x=solve_ivp(motored,ts,x0,method='DOP853').y
        
        P1[l+1]=x[0][-1]
        T1[l+1]=x[1][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
    

    x0 = [P0,T0,Qa0,Qp0,W0]
    #fechamento da valvula ao inicio da combustao
    for i in range(len(ti)-1): #range(0,98)
        ts=[ti[i],ti[i+1]]
        x=solve_ivp(motor,ts,x0,args=(0,),method='DOP853').y
        ke = KReag(T[i])
        k_eval.append(ke)
        P[i+1]=x[0][-1] #ultima linha,coluna1
        T[i+1]=x[1][-1]
        Qa[i+1]=x[2][-1]
        Qp[i+1]=x[3][-1]
        W[i+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        massFraction[i+1] = xt

    #combustao    
    for j in range (len(ti)-1,len(ti)+len(tj)-1):
        ts=[tk[j],tk[j+1]]
        
        xt = wiebe(tk[j])
        dxdt = dx(tk[j])
        
        ke = KReag(T[j])*(1-xt)+KProd(T[j])*xt
        k_eval.append(ke)
        
        x=solve_ivp(motorComb,ts,x0,args=(P1[j],xt,dxdt,0),method='DOP853').y

        P[j+1]=x[0][-1]
        T[j+1]=x[1][-1]
        Qa[j+1]=x[2][-1]
        Qp[j+1]=x[3][-1]
        W[j+1]=x[4][-1]
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        massFraction[j+1] = xt
    
    for f in range(len(tj)+len(ti)-1,len(tk)-1): #200 ao 298
        ts=[tk[f],tk[f+1]]
        
        xt = wiebe(tk[f])
        dxdt = dx(tk[f])
        #dxdt = 0
        
        ke = KReag(T[f])*(1-xt)+KProd(T[f])*xt
        k_eval.append(ke)
        
        x=solve_ivp(motorComb,ts,x0,args=(P1[f],xt,dxdt,0),method='DOP853').y
        
        P[f+1]=x[0][-1]
        T[f+1]=x[1][-1]
        Qa[f+1]=x[2][-1]
        Qp[f+1]=x[3][-1]
        W[f+1]=x[4][-1]
        
        
        x0=[x[0][-1],x[1][-1],x[2][-1],x[3][-1],x[4][-1]]
        
        massFraction[f+1] = xt
        
    return P,T,Qa,Qp,W,massFraction,k_eval