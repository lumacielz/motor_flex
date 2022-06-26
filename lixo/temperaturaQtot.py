#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jan  1 16:13:24 2022

@author: luiza.maciel
"""


from scipy.optimize import minimize
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from decimal import  Decimal
from math import radians
import numpy as np
from sympy.core.numbers import comp
import sys

file2='./tabelaHS.xlsx'
df2=pd.read_excel(file2)

T0=float(input("Temperatura dos reagentes: "))
P0=float(input("Pressão[atm]:  "))
#estimateT=float(input("Estimativa da temperatura final: "))
combustivel=input("combustível [gasolina,etanol,GNV,H2,outro]")
Ts=sp.symbols("Ts")

def stringToFloat(inputLista):
    lista=[]
    for i in inputLista.split(","):
        lista.append(float(i))
    return lista

def composicaoCombustivel(c):
    coef=c.split("C")[1].split("H")
    coef[1:]=(coef[1].split("O"))
    return list(map(lambda x: float(x),coef))

r=11
L=0.144
D=0.081
S=0.086

R=S/2 #raio do virabrequim
teta=sp.symbols("teta")
s=R*sp.cos(teta)+sp.sqrt(L**2-R**2*(sp.sin(teta))**2)  #deslocamento

Vo=(3.14*D**2/4)*(L+R-s+2*R/(r-1)) #volume 
print(Vo.subs(teta,0))

teta0d=-164
teta_combd=-4
delta_tetad=38
tetafd=146
teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)

V0=Vo.subs(teta,teta0)
V=V0

while True:
    if combustivel=="gasolina":
    #combustivel equivalente
    #    carb=6.67
    #    h=12.8
    #    o=0.533  
    #    n=0
    #    OC=9.6
    #    complete_combustion_H2O=6.4
    #    complete_combustion_CO2=6.67
        
        c,h,o,n=8*0.5462+2*0.4538,18*0.5462+0.4538*6,0.4538,0
        
        a0=-8102.9675849653
        a1=7963.204888950
        a2=-2923.24238668569
        a3=509.144026930886
        a4=-42.2764373650608
        a5=1.35175803975684
        Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
        
        h0_comb=0.5462*(-208450)+0.4538*(-235310)
        s0_comb=0.5462*(466.514)+0.4538*(282.444)
        mcomb=(12*c+h+16*o)
        
        PCI=39249*10**3
        
        break
    
    elif combustivel=="alcool":
        #combustivel equivalente
    #    OC=3.2
    #    carb=2.15
    #    h=6.62
    #    o=1.23
    #    complete_combustion_CO2=2.15
    #    complete_combustion_H2O=3.31
        
        c,h,o,n=1.709,5.418801,0
    
        a0=-12482.8740213179
        a1=10263.4623453335
        a2=-3316.7850402598
        a3=526.309291795851
        a4=-40.8869367350809
        a5=1.24555084441151
        Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
        
        h0_comb=0.8547*(-235310)+0.1453*(-241826)
        s0_comb=0.8547*(282.444)+0.1453*(188.835)
        mcomb=(12*c+h+16*o)
        
        PCI=24804*10**3
        
        break
    
    elif combustivel=="GNV":
        met=0.92285
        et=5.455*10**-2
        prop=1.115*10**-2
        but=0.125*10**-2
        co2=0.255*10**-2
        r=0.765*10**-2
        c=met+2*et+3*prop+4*but+co2
        h=met*4+et*6+prop*8+but*10
        o=co2*2
        n=r
    
        a0=-12713.9307245771
        a1=10272.7734235011
        a2=-3251.2253041433
        a3=504.476942937525
        a4=-38.3629908265519
        a5=1.14655193729902
        Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
        
        h0_met=-74873
        s0_met=186.251
        h0_et=-84740
        s0_et=229.597
        h0_prop=-103900
        s0_prop=269.917
        h0_but=-126200
        s0_but=306.647
        h0_co2=-393522
        s0_co2=213.795
        h0_n2=0
        s0_n2=191.609
        h0_comb=met*h0_met+et*h0_et+prop*h0_prop+but*h0_but+co2*h0_co2+r*h0_n2
        s0_comb=met*s0_met+et*s0_et+prop*s0_prop+but*s0_but+co2*s0_co2+r*s0_n2
        mcomb=(met*16+et*30+prop*44+but*58+co2*44+28*r)
        
        PCI=48737*10**3
        
        break 
    
    elif combustivel == "H2":
        c,h,o,n=0,2,0,0
        mcomb=2
        
        PCI=28700*4186.8
        
        break
    
    elif combustivel == "outro":
        molecularComposition = input("digite a composição do combustível na forma CnHnOn: ")
        c,h,o,n=composicaoCombustivel(molecularComposition),0
        mcomb=c*12+h+o*16
        h0_comb=float(input("entalpia de formação do combustivel: "))
        s0_comb=float(input("entropia de formação do combustivel: "))
        a0,a1,a2,a3,a4,a5=stringToFloat(input("coeficientes da equação logaritimica de Cp em funcao da temperatura: "))
        Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5  #J/molK
        
        PCI=float(input("poder calorifico do combustivel: "))
        
    else:
        print("digite um combustível válido!!!")


mols=float(input("quantidade de combustivel: "))
oxid=input("[0] ar atmosferico [1] somente O2")
OC=float(input("quantidade de oxigenio: "))
NI=0
if oxid == "0":
    NI=3.76*OC
estimative=stringToFloat(input("estimativa inicial para a composicao de equilibrio [Combustivel,O2,H2O,CO2,CO,O,N,NO2,NO,OH,H2,H,N2]: "))

#comb,o2,h2o,co2,co,o,n,no2,no,oh,h2,h,n2
m=[mcomb,32,18,44,28,16,14,46,30,17,2,1,28]

total_mass=mols*mcomb+32*OC+NI*28

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
   
    Ntot=sum(composicoes)
    
    deltaS=interpolacaoS(T)
    constant=Ru*T
    
    h_comb,h_o2, h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2=interpolacaoH(T)
    
    gt_comb=(h_comb-T*deltaS["deltasComb"])/constant
    gt_o2=(h_o2-T*deltaS["deltasO2"])/constant
    gt_h2o=(h_h2o-T*deltaS["deltasH2O"])/constant
    gt_co2=(h_co2-T*deltaS["deltasCO2"])/constant
    gt_co=(h_co-T*deltaS["deltasCO"])/constant
    gt_o = (h_o - T * deltaS["deltasO"]) / constant
    gt_n = (h_n - T * deltaS["deltasN"]) / constant
    gt_no2 = (h_no2 - T * deltaS["deltasNO2"]) / constant
    gt_no = (h_no - T * deltaS["deltasNO"]) / constant
    gt_oh = (h_oh - T * deltaS["deltasOH"]) / constant
    gt_h2 = (h_h2 - T * deltaS["deltasH2"])/ constant
    gt_h = (h_h - T * deltaS["deltasH"]) / constant
    gt_n2=(h_n2-T*deltaS["deltasN2"])/constant

    gt=[gt_comb,gt_o2,gt_h2o,gt_co2,gt_co,gt_o,gt_n,gt_no2,gt_no,gt_oh,gt_h2,gt_h,gt_n2]
    deltaG=0
    for n,g in zip(composicoes,gt):
        if n!=0:
            deltaG+=n * (g + float((Decimal(n) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))
   
    return (deltaG)

def restricaoC(composicoes):
    return mols*c-composicoes[0]*c-composicoes[3]-composicoes[4]
def restricaoH(composicoes):
    return mols*h-composicoes[0]*h-composicoes[2]*2-composicoes[9]-composicoes[10]*2-composicoes[11]
def restricaoO(composicoes):
    return mols*o+2*OC-composicoes[0]*o-2*composicoes[1]-composicoes[2]-composicoes[3]*2-composicoes[4]-composicoes[5]-composicoes[7]*2-composicoes[8]-composicoes[9]
def restricaoN(composicoes):
    return mols*n+2*NI-n*composicoes[0]-composicoes[12]*2-composicoes[6]-composicoes[7]-composicoes[8]

def interpolacaoH(i):
    for k in range(36):
        if i>=df2['T'][k] and i<=df2['T'][k+1]:
            deltahH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k]
            deltahCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k]
            deltahO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k]
            deltahCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k]
            deltahN2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hN2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hN2'][k]
            deltahNO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hNO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hNO'][k]
            deltahNO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hNO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hNO2'][k]
            deltahOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k]
            deltahH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k]
            deltahH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k]
            deltahN=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hN'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hN'][k]
            deltahO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k]
            
    
    h_o2 = h0_o2 + deltahO2    
    h_h2o = h0_h2o + deltahH2O
    h_co2 = h0_co2 + deltahCO2
    h_co = h0_co + deltahCO
    h_o = h0_o + deltahO
    h_n = h0_n + deltahN
    h_no2 = h0_no2 + deltahNO2
    h_no = h0_no + deltahNO
    h_oh = h0_oh + deltahOH
    h_h2 = h0_h2 + deltahH2
    h_h = h0_h + deltahH
    h_n2 = h0_n2 + deltahN2

    if combustivel=="H2":
        h_comb=h_h2
        h_h2=0
    else:
        h_comb = float(h0_comb + sp.integrate(Cp, (Ts, 298.15, i)))
    
    return h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2

def interpolacaoS(i):
    for k in range(36):
            if i>=df2['T'][k] and i<=df2['T'][k+1]:
                deltasH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k]
                deltasCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k]
                deltasO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k]
                deltasCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO'][k]
                deltasN2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sN2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sN2'][k]
                deltasNO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sNO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sNO'][k]
                deltasNO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sNO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sNO2'][k]
                deltasOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sOH'][k]
                deltasH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2'][k]
                deltasH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH'][k]
                deltasN=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sN'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sN'][k]
                deltasO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO'][k]
                if combustivel=="H2":
                    deltasComb=deltasH2
                    deltasH2=0
                else:
                    deltasComb=float(s0_comb + sp.integrate(Cp/Ts,(Ts,298.15,i)))
                
    return {"deltasComb":deltasComb,"deltasO2":deltasO2,"deltasH2O":deltasH2O,"deltasCO2":deltasCO2,"deltasCO":deltasCO,"deltasO":deltasO,"deltasN":deltasN,"deltasNO2":deltasNO2,"deltasNO":deltasNO,"deltasOH":deltasOH,"deltasH2":deltasH2,"deltasH":deltasH,"deltasN2":deltasN2}

def proxTemperatura(composicao,T0,T1):
    
    Ncomb=composicao[0]
    NO2=composicao[1]
    NH2O=composicao[2]
    NCO2=composicao[3]
    NCO=composicao[4]
    NO=composicao[5]
    NN=composicao[6]
    NNO2=composicao[7]
    NNO=composicao[8]
    NOH=composicao[9]
    NH2=composicao[10]
    NH=composicao[11]
    NN2=composicao[12]
    
    while True:
        h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2=interpolacaoH(T0)
        entalpy_products0=Ncomb*h_comb+NO2*h_o2+NH2O*h_h2o+NCO2*h_co2+NCO*h_co+NO*h_o+NN*h_n+NNO2*h_no2+NNO*h_no+NOH*h_oh+NH2*h_h2+NH*h_h+NN2*h_n2
        zero0=entalpy_products0-entalpy_reactants+Qtot
        
        h_comb,h_o2,h_h2o,h_co2,h_co,h_o,h_n,h_no2,h_no,h_oh,h_h2,h_h,h_n2=interpolacaoH(T1)
        entalpy_products1=Ncomb*h_comb+NO2*h_o2+NH2O*h_h2o+NCO2*h_co2+NCO*h_co+NO*h_o+NN*h_n+NNO2*h_no2+NNO*h_no+NOH*h_oh+NH2*h_h2+NH*h_h+NN2*h_n2
        zero1=entalpy_products1-entalpy_reactants+Qtot
        
        print(zero1,zero0)
        nextT=T1-zero1*(T1-T0)/(zero1-zero0)
        
        if abs(zero1)<0.2:
            break
        else:
            T0=T1
            T1=nextT
    return nextT


entalpy_reactants=interpolacaoH(T0)[0]*mols+OC*interpolacaoH(T0)[1]+NI*interpolacaoH(T0)[12]

cons=[{'type': 'eq', 'fun': restricaoC},{'type': 'eq', 'fun': restricaoO},{'type': 'eq', 'fun':restricaoH},{'type': 'eq', 'fun': restricaoN}]
#estimative=[10**-10,10**-10,2,0,0,10**-10,0,0,0,10**-10,0,10**-10,0]
#estimative=[1e-10,1e-10,2,2,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10,1e-10]

def bounds (composicao):
    bnds=[]
    for c in composicao:
        if c==0:
            bnds.append((0,sys.float_info.min))
        else:
            bnds.append((0,None))
    return tuple(bnds)
bnds=bounds(estimative)

P=P0
T=T0
mi=mols*m[0]
Qt=0.87*mi*PCI

for t in np.linspace(teta_comb,teta_comb+delta_teta,100):
    V=Vo.subs(teta,t)
    
    sol = minimize(gibbs, estimative,args=(T,P), method='SLSQP', bounds=bnds, constraints=cons)
    print(sol)
    total_mols=sum(sol.x)
    mc=sol.x[0]*m[0]
    xt=(mi-mc)/mi
    print(xt)
    Qtot=Qt*xt*10**-3
    
    nextT=proxTemperatura(sol.x,T,T+10)
    print(nextT,xt)
    
    T=nextT
    P=total_mols*Ru*T/V*0.00986923
    
    
    
    
total_mols=sum(sol.x)


print("COMPOSICAO EQUILIBRIO: " +  str({'Combustivel': sol.x[0] , 'O2': sol.x[1],
       'H2O': sol.x[2] , 'CO2': sol.x[3] , 'CO': sol.x[4] ,'O':sol.x[5] ,'N':sol.x[6] ,'NO2':sol.x[7] ,'NO':sol.x[8],'OH':sol.x[9] ,'H2':sol.x[10] ,'H':sol.x[11] ,'N2':sol.x[12]}))
print({'MASS_FRACTION':{'Combustivel': sol.x[0] * m[0] / total_mass, 'O2': sol.x[1] * m[1] / total_mass,
       'H2O': sol.x[2] * m[2] / total_mass, 'CO2': sol.x[3] * m[3] / total_mass, 'CO': sol.x[4] * m[4] / total_mass,'O':sol.x[5] * m[5] / total_mass,'N':sol.x[6] * m[6] / total_mass,'NO2':sol.x[7] * m[7] / total_mass,'NO':sol.x[8] * m[8] / total_mass,'OH':sol.x[9] * m[9] / total_mass,'H2':sol.x[10] * m[10] / total_mass,'H':sol.x[11] * m[11] / total_mass,'N2':sol.x[12] * m[12] / total_mass}})
print({'MOLAR_FRACTION':{'Combustivel': sol.x[0] / total_mols, 'O2': sol.x[1] / total_mols, 'H2O': sol.x[2] / total_mols, 'CO2':sol.x[3]/total_mols, 'CO':sol.x[4]/total_mols,'O':sol.x[5]/ total_mols,'N':sol.x[6]/total_mols,'NO2':sol.x[7]/total_mols,'NO':sol.x[8]/ total_mols,'OH':sol.x[9]/ total_mols,'H2':sol.x[10]/ total_mols,'H':sol.x[11]/ total_mols,'N2':sol.x[12]/ total_mols}})

