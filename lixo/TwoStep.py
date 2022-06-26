#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov  6 12:07:43 2021

@author: luiza.maciel
"""
from math import exp,radians

import sympy as sp
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint,solve_ivp

teta=sp.symbols("teta")
file1='./dados.txt'

arq=open(file1)

data=arq.readlines()
arq.close()
Temp=[]
Press=[]
for t in data[3:303]:
    Temp.append(float(t))
Pressure=[]
for p in data [304:]:
    Press.append(float(p))

A=(5.7*10**11)
E=(15098)
m=0.25
n=1.5
o=8.5

O0=(o/((1+o)+3.76*o))*(Press[100]*10**-3/(8315*Temp[100]))
C0=(1/((1+o)+3.76*o))*(Press[100]*10**-3/(8315*Temp[100]))
Concentrations0=[C0,12.5*C0,0,0,0]
C1=C0
m_c=8*12+18
m_t=m_c+o*32+3.76*o*28
m0=1/((1+o)+3.76*o)*m_c

teta0d= -164
tetafd=146
teta_combd=-4
delta_tetad=39
teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)

ti=np.linspace(teta0,teta_comb,100) #fechamento da admissao ate inicio da combustao
tj=np.linspace(teta_comb,teta_comb+delta_teta,100) #combustao
tf=np.linspace(teta_comb+delta_teta,tetaf,100) #fim da combustao ate abertura do escapamento

tk=np.ones(len(ti)+len(tj)+len(tf)) 
tk[:len(ti)]=ti #0 a 99
tk[len(ti):len(tj)+len(ti)]=tj #100 ao 199

tk[len(tj)+len(ti):len(tf)+len(tj)+len(ti)]=tf #200 ao 299
N=2500*6.28/60

def dc (Concentrations,t,T):
    C=Concentrations[0]
    O2=Concentrations[1]
    CO=Concentrations[2]
    CO2=Concentrations[3]
    H2O=Concentrations[4]
    
    ka=A*exp(-15098/T)*np.sign(C)*abs(C)**m*O2**n
    kbf=10**14.6*exp(-20131/T)*CO*H2O**0.5*O2**0.25
    kbr=5*10**8*exp(-20131/T)*CO2
    
    dC=-ka/N
    dO2=(-o*ka-0.5*8*kbf+8*0.5*kbr)/N
    dCO=(-kbf+ka+kbr)*8/N
    dCO2=(kbf-kbr)*8/N
    dH2O=9*ka/N
    return dC,dO2,dCO,dCO2,dH2O

lista=[0]

for t in range(len(tk))[100:199]:
    T=Temp[t]
    x=odeint(dc,Concentrations0,[tk[t],tk[t+1]],args=(T,),atol=10**-17,rtol=10**-17)
    Concentrations0=x[-1]
    lista.append((m0-(x[-1][0]*8315*Temp[t]/(Press[t]*10**-3))*m_c)/m0)
    print(x[-1][0])

time=np.linspace(teta_combd,teta_combd+delta_tetad,100)
plt.scatter(time,lista)