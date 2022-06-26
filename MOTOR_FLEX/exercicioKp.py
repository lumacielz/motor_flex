#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 13:46:20 2022´

@author: luiza.maciel
"""

from math import exp
#Kp da reação 2CO2 -> 2CO + O2
R=8.314
#T = 500
#gco =(-110527+5932) - T * 212.833
#gco2=(-393522+8305) - T * 234.902
#go2 = 6086 - T * 220.693
#deltaG = 2*gco + go2 - 2*gco2
#
#lnKp = -deltaG/(R*T)
#
#
##CH4 + 2O2 -> CO2+2H2O
#T=298.15
#gch4 = -74873 - T*186.251
#go2 = -T*205.148
#gco2 = -393522 - T*213.794
#gh2o = -241826 - T*188.835
#
#deltaG = gco2+2*gh2o - gch4-2*go2
#lnKp = -deltaG/(R*T)
#Kp = exp(lnKp)

import pandas as pd
from math import exp

file='./tabelaHS.xlsx'
df2=pd.read_excel(file)
#entalpias de formacao 
h0_o2=0
h0_h2o=-241820
h0_co2=-393520
h0_co=-110530
h0_o = 249190
h0_n = 472650
h0_no2 = 33100
h0_no = 90291
h0_oh = 39460
h0_h2 = 0
h0_h = 218000
h0_n2=0

def interpolacaoH(i):
    for k in range(36):
        if i>=df2['T'][k] and i<=df2['T'][k+1]:
            deltahCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k]
            deltahO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k]
            deltahCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k]
            
    h_o2 = h0_o2 + deltahO2    
    h_co2 = h0_co2 + deltahCO2
    h_co = h0_co + deltahCO
    
    return {"O2":h_o2,"CO2":h_co2,"CO":h_co}

def interpolacaoS(i):
    for k in range(36):
            if i>=df2['T'][k] and i<=df2['T'][k+1]:
                deltasCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k]
                deltasO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k]
                deltasCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO'][k]
                
    return {"O2":deltasO2,"CO2":deltasCO2,"CO":deltasCO}

def Kpf(t):
    h = interpolacaoH(t)
    ds = interpolacaoS(t)
    gco = h["CO"] - t * ds["CO"]
    gco2 = h["CO2"] - t * ds["CO2"]
    go2 = h["O2"] - t * ds["O2"]
    deltaG = gco2 - gco - 0.5*go2
    
    Kp = exp(-deltaG/(R*t))
    Kc = Kp/((R*9.86923e-6*t))**-0.5
    return Kp


