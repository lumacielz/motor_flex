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

A=(4.6*10**11)
E=(15098)
m=0.25
n=1.5
o=12.5

O0=(o/((1+o)+3.76*o))*(Press[100]*10**-3/(8315*Temp[100]))
C0=(1/((1+o)+3.76*o))*(Press[100]*10**-3/(8315*Temp[100]))
print(C0)
C1=C0
m_c=8*12+18
m_t=m_c+o*32+3.76*o*28
teta0d= -164
tetafd=157
teta_combd=-4
delta_tetad=37
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

def dc (C,t,T):
    #O2=O0+(C-C0)/12.5
    O2=12.5*C
    dC=-(A*exp(-E/T)*np.sign(C)*abs(C)**m*(O2)**n)/N
    return dC
m0=(1/((1+o)+3.76*o))*m_c
lista=[(m-m)/m]
print(lista)
for t in range(len(tk))[100:199]:
    T=Temp[t]
    x=odeint(dc,C1,[tk[t],tk[t+1]],args=(T,),rtol=1e-18,atol=1e-18)
    C1=x[1]
    lista.append((m0-(x[-1][0]*8315*Temp[t]/(Press[t]*10**-3))*m_c)/m0)
    print(x[1])

time=np.linspace(teta_combd,teta_combd+delta_tetad,100)
plt.plot(time,lista)
plt.ylabel("Mass Fraction Burned")