#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
"""


from scipy.optimize import minimize
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from decimal import  Decimal
from math import radians
import numpy as np

file2='./tabelaHS.xlsx'
df2=pd.read_excel(file2)

T0=298.15
P0=22
estimateT=3000

def composicaoCombustivel(c):
#    coef=[0,0,0,0]
#    if c.find("C") !=-1:
#        coef[0]=c.split("C")[1]
#    if c.find("H") != -1:
#        coef[1]=c.split("H")[1]
#    if c.find("O") != -1:
#        coef[2]=c.split("O")[1]
#    if c.find("N") != -1:
#        coef[3]=c.split("N")[1]
        
    coef=c.split("C")[1].split("H")
    coef[1:]=(coef[1].split("O"))
    return list(map(lambda x: float(x),coef))

#comb=input("digite a composição do combustível na forma CnHnOn: ")
comb=composicaoCombustivel("C0H2O0")
mols=float(input("quantidade de combustivel: "))
oxid=float(input("quantidade de oxigenio: "))
#se for nulo sera o esteq

c=comb[0]
h=comb[1]
o=comb[2]
mcomb=c*12+h+o*16

m=[mcomb,32,18,17,16,1]

total_mass=mols*mcomb+32*oxid

Ru=8.314

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

def gibbs (composicoes,T,P):
    
    NH2=composicoes[0]
    NO2=composicoes[1]
    NH2O=composicoes[2]
    NOH=composicoes[3]
    NO=composicoes[4]
    NH=composicoes[5]
    Ntot=sum(composicoes)
    
    deltaS=interpolacaoS(T)
    constant=Ru*T
    
    h_h2,h_o2,h_h2o,h_oh,h_o,h_h=interpolacaoH(T)

    gt_o2=(h_o2-T *deltaS["deltasO2"])/constant
    gt_h2o=(h_h2o-T*deltaS["deltasH2O"])/constant
    gt_oh = (h_oh - T* deltaS["deltasOH"]) / constant
    gt_h2 = (h_h2 - T * deltaS["deltasH2"])/ constant
    gt_h = (h_h - T * deltaS["deltasH"]) / constant
    gt_o = (h_o - T * deltaS["deltasO"]) / constant
   
    deltaG=NH2 * (gt_h2 + float((Decimal(NH2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) +\
       NO2 * (gt_o2 + float((Decimal(NO2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))+\
       NH2O * (gt_h2o + float((Decimal(NH2O) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) +\
       NOH * (gt_oh + float((Decimal(NOH) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NO * (gt_o + float((Decimal(NO) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NH * (gt_h + float((Decimal(NH) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))
   
    return (deltaG)

def restricaoO(composicoes):
    return (2*oxid)-2*composicoes[1]-composicoes[2]-composicoes[3]-composicoes[4]
def restricaoH(composicoes):
    return h*mols-2*composicoes[0]-2*composicoes[2]-composicoes[3]-composicoes[5]

def interpolacaoH(i):
    for k in range(36):
            if i>=df2['T'][k] and i<=df2['T'][k+1]:
                deltahH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k]
                deltahO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k]
                deltahH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k]
                deltahOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k]
                deltahO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k]
                deltahH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k]
                
    deltaH0={"deltahH2":deltahH2,"deltahO2":deltahO2,"deltahH2O":deltahH2O,"deltahOH":deltahOH,"deltahO":deltahO,"deltahH":deltahH}
        
    h_h2o = h0_h2o + deltaH0["deltahH2O"]
    h_o2 = h0_o2 + deltaH0["deltahO2"]
    h_oh = h0_oh + deltaH0["deltahOH"]
    h_h2 = h0_h2 + deltaH0["deltahH2"]
    h_h = h0_h + deltaH0["deltahH"]
    h_o = h0_o + deltaH0["deltahO"]
    
    return h_h2,h_o2,h_h2o,h_oh,h_o,h_h

def interpolacaoS(i):
    for k in range(36):
            if i>=df2['T'][k] and i<=df2['T'][k+1]:
                deltasH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2'][k]
                deltasH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k]
                deltasO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k]
                deltasOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sOH'][k]
                deltasO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO'][k]
                deltasH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH'][k]
          
                
    return {"deltasH2":deltasH2,"deltasO2":deltasO2,"deltasH2O":deltasH2O,"deltasOH":deltasOH,"deltasO":deltasO,"deltasH":deltasH}
def temperaturaAdiabatica(composicao,T0,T1):
    NH2=composicao[0]
    NO2=composicao[1]
    NH2O=composicao[2]
    NOH=composicao[3]
    NO=composicao[4]
    NH=composicao[5]
   
    while True:
        h_h2,h_o2,h_h2o,h_oh,h_o,h_h=interpolacaoH(T0)
        entalpy_products0=NH2*h_h2+NO2*h_o2+NH2O*h_h2o+NOH*h_oh+NO*h_o+NH*h_h
        zero0=entalpy_products0-entalpy_reactants
        
        h_h2,h_o2,h_h2o,h_oh,h_o,h_h=interpolacaoH(T1)
        entalpy_products1=NH2*h_h2+NO2*h_o2+NH2O*h_h2o+NOH*h_oh+NO*h_o+NH*h_h
        zero1=entalpy_products1-entalpy_reactants
        
        nextT=T1-zero1*(T1-T0)/(zero1-zero0)
        print(nextT)
        if abs(zero1)<0.5:
            break
        else:
            T0=T1
            T1=nextT
    return nextT

def composicaoAdiabatica(estimative,T0,P0):
    P=P0
    T=T0
    while True:
        
        sol = minimize(gibbs, estimative,args=(T,P), method='SLSQP', bounds=bnds, constraints=cons)
        
        total_mols = sum(i for i in sol.x)
        
        nextTemp=float(temperaturaAdiabatica(sol.x,T,T+2))
        
        if abs(nextTemp-T)<0.2:
            break
        else:
            T=nextTemp
            estimative=sol.x
            
    return nextTemp,P0,sol.x,total_mols

entalpy_reactants=interpolacaoH(T0)[0]*mols+oxid*interpolacaoH(T0)[1]

bnds=((0,None),(0,None),(0,None),(0,None),(0,None),(0,None))
cons=[{'type': 'eq', 'fun': restricaoO},{'type': 'eq', 'fun':restricaoH}]
estimative=[10**-10,10**-10,2,10**-10,10**-10,10**-10]
P=P0
T=estimateT
sol = minimize(gibbs, estimative,args=(T,P), method='SLSQP', bounds=bnds, constraints=cons)

print({"H2":sol.x[0],"O2":sol.x[1],"H2O":sol.x[2],"OH":sol.x[3],"O":sol.x[4],"H":sol.x[5]})
#temperaturaAdiabatica(sol.x,T,T+2)
adiab=composicaoAdiabatica(sol.x,T,P)
print({"Tad":adiab[0],"composicaoAdiab":{"H2":adiab[2][0],"O2":adiab[2][1],"H2O":adiab[2][2],"OH":adiab[2][3],"O":adiab[2][4],"H":adiab[2][5]}})