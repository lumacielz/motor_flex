#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 13 08:20:11 2021

@author: luiza.maciel
"""
#Proxima temperatura e atualizada para Tadiabatica

import sys
from scipy.optimize import minimize
import sympy as sp
import pandas as pd
import matplotlib.pyplot as plt
from decimal import  Decimal
from math import radians
import numpy as np

file1='./dados.txt'
file2='./tabelaHS.xlsx'
df2=pd.read_excel(file2)
arq=open(file1)

data=arq.readlines()
arq.close()
#pegando os dados  de temperatura e pressão 
m_c=float(data[1])
T0=float(input("Temperatura inicial: "))
P0=float(input("Pressão inicial: "))
#estimative=input("Estimativa inicial [Combustivel,O2,H2O,CO2,CO,O,N,NO2,NO,OH,H2,H,N2]: ")
#estimative0=[float(x) for x in estimative.split(",")]
T0=300
P0=67500
Temp=[T0]
Press=[P0]

teta=sp.symbols("teta")
L=0.144
D=0.081
R=D/2
r=11
s=R*sp.cos(teta)+sp.sqrt(L**2-R**2*(sp.sin(teta))**2)  #deslocamento

Vo=(3.14*D**2/4)*(L+R-s+2*R/(r-1)) #volume 

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

V0=Vo.subs(teta,teta0)
combustivel=input('COMBUSTIVEL: ')
#definindo composição , Cp's e  entalpias de formação dos combustíveis

if combustivel=="gasolina":
    
#    carb=8*0.5462+2*0.4538
#    h=18*0.5462+0.4538*6
#    o=0.4538
#    n=0
    carb=8
    h=18
    o=0
    n=0
    
    a0=-8102.9675849653
    a1=7963.204888950
    a2=-2923.24238668569
    a3=509.144026930886
    a4=-42.2764373650608
    a5=1.35175803975684
    
    h0_comb=-208450
    s0_comb=466.514
#    h0_comb=0.5462*(-208450)+0.4538*(-235310)
#    s0_comb=0.5462*(466.514)+0.4538*(282.444)
    mcomb=(12*carb+h+16*o)

elif combustivel=="alcool":

    carb=1.709
    o=1
    h=5.4188
    n=0

    a0=-12482.8740213179
    a1=10263.4623453335
    a2=-3316.7850402598
    a3=526.309291795851
    a4=-40.8869367350809
    a5=1.24555084441151
    
    h0_comb=0.8547*(-235310)+0.1453*(-241826)
    s0_comb=0.8547*(282.444)+0.1453*(188.835)
    mcomb=(12*carb+h+16*o)

elif combustivel=="GNV":
    met=0.92285
    et=5.455*10**-2
    prop=1.115*10**-2
    but=0.125*10**-2
    co2=0.255*10**-2
    r=0.765*10**-2
    carb=met+2*et+3*prop+4*but+co2
    h=met*4+et*6+prop*8+but*10
    o=co2*2
    n=r

    a0=-12713.9307245771
    a1=10272.7734235011
    a2=-3251.2253041433
    a3=504.476942937525
    a4=-38.3629908265519
    a5=1.14655193729902
    
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
    
elif combustivel=="H2":
    h=2
    carb=0
    o=0
    mcomb=2
    h0_comb=0
    
else:
    print("digite um combustível válido!!!")

#calcula coeficientes estequiometricos
complete_combustion_CO2=carb
complete_combustion_H2O=h/2
OC=(complete_combustion_H2O+2*complete_combustion_CO2-o)/2

#massas moleculares de cada composto e eq do Cp
T=sp.symbols('T')
m=[mcomb,32,18,44,28,16,14,46,30,17,2,1,28]
mols=m_c/m[0]
total_mass=mols*(1*m[0]+OC*m[1]+OC*3.76*m[12])

Cp=a0+a1*sp.log(T)+a2*sp.log(T)**2+a3*sp.log(T)**3+a4*sp.log(T)**4+a5*sp.log(T)**5  #J/molK
if combustivel=="H2":
    Cp=0
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

#armazenamento da evolução da reacao
composicao_comb=[]
composicao_o2=[]
composicao_h2o=[]
composicao_co2=[]
composicao_co=[]
composicao_o=[]
composicao_n=[]
composicao_no2=[]
composicao_no=[]
composicao_oh=[]
composicao_h2=[]
composicao_h=[]
composicao_n2=[]

composicao_comb.append(mols * m[0] / total_mass)
composicao_o2.append(mols*OC * m[1] / total_mass)
composicao_h2o.append(0 * m[2] / total_mass)
composicao_co2.append(0 * m[3] / total_mass)
composicao_co.append(0 * m[4] / total_mass)
composicao_o.append(0 * m[5] / total_mass)
composicao_n.append(0* m[6] / total_mass)
composicao_no2.append(0 * m[7] / total_mass)
composicao_no.append(0 * m[8] / total_mass)
composicao_oh.append(0 * m[9] / total_mass)
composicao_h2.append(0 * m[10] / total_mass)
composicao_h.append(0* m[11] / total_mass)
composicao_n2.append(mols*3.76*OC * m[12] / total_mass)


def gibbs (composicoes,T,P):
    
    deltaS=interpolacaoS(T)
    constant=Ru*T
    
    h_h2o,h_co2,h_o2,h_co,h_n2,h_comb,h_no,h_no2,h_oh,h_h2,h_h,h_n,h_o=interpolacaoH(T)

    gt_comb=(h_comb-T*deltaS["deltasComb"])/constant
    gt_o2=(h_o2-T*deltaS["deltasO2"])/constant
    gt_h2o=(h_h2o-T*deltaS["deltasH2O"])/constant
    gt_co2=(h_co2-T*deltaS["deltasCO2"])/constant
    gt_co=(h_co-T*deltaS["deltasCO"])/constant
    gt_n2=(h_n2-T*deltaS["deltasN2"])/constant
    gt_no = (h_no - T * deltaS["deltasNO"]) / constant
    gt_no2 = (h_no2 - T * deltaS["deltasNO2"]) / constant
    gt_oh = (h_oh - T * deltaS["deltasOH"]) / constant
    gt_h2 = (h_h2 - T * deltaS["deltasH2"])/ constant
    gt_h = (h_h - T * deltaS["deltasH"]) / constant
    gt_n = (h_n - T * deltaS["deltasN"]) / constant
    gt_o = (h_o - T * deltaS["deltasO"]) / constant
    
    Ncomb=composicoes[0]
    NO2=composicoes[1]
    NH2O=composicoes[2]
    NCO2=composicoes[3]
    NCO=composicoes[4]
    NO=composicoes[5]
    NN=composicoes[6]
    NNO2=composicoes[7]
    NNO=composicoes[8]
    NOH=composicoes[9]
    NH2=composicoes[10]
    NH=composicoes[11]
    NN2=composicoes[12]
    Ntot=sum(composicoes)
   
    deltaG=Ncomb * (gt_comb + float((Decimal(Ncomb) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NCO2 * (gt_co2 + float((Decimal(NCO2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NH2O * (gt_h2o + float((Decimal(NH2O) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) +\
       NO2 * (gt_o2 + float((Decimal(NO2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))+\
       NN2 * (gt_n2 + float((Decimal(NN2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))+\
       NCO * (gt_co + float((Decimal(NCO) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NO * (gt_o + float((Decimal(NO) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NN * (gt_n + float((Decimal(NN) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NNO2 * (gt_no2 + float((Decimal(NNO2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NNO * (gt_no + float((Decimal(NNO) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NOH * (gt_oh + float((Decimal(NOH) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NH2 * (gt_h2 + float((Decimal(NH2) * Decimal("%.15f" % P) / Decimal(Ntot)).ln())) + \
       NH * (gt_h + float((Decimal(NH) * Decimal("%.15f" % P) / Decimal(Ntot)).ln()))
   
    return (deltaG)

def restricaoC(composicoes):
    return mols*carb-composicoes[0]*carb-composicoes[3]-composicoes[4]
def restricaoO(composicoes):
    return mols*(o+2*OC)-composicoes[0]*o-2*composicoes[1]-composicoes[2]-composicoes[4]-composicoes[3]*2-composicoes[5]-2*composicoes[7]-composicoes[8]-composicoes[9]
def restricaoH(composicoes):
    return mols*h-h*composicoes[0]-2*composicoes[2]-composicoes[9]-2*composicoes[10]-composicoes[11]
def restricaoN(composicoes):
    return mols*(OC*3.76*2+n)-n*composicoes[0]-composicoes[12]*2-composicoes[6]-composicoes[7]-composicoes[8]

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
            
    deltaH0={"deltahH2O":deltahH2O,"deltahCO2":deltahCO2,"deltahO2":deltahO2,"deltahCO":deltahCO,"deltahN2":deltahN2,"deltahNO":deltahNO,"deltahNO2":deltahNO2,"deltahOH":deltahOH,"deltahH2":deltahH2,"deltahH":deltahH,"deltahN":deltahN,"deltahO":deltahO}
        
    h_h2o = h0_h2o + deltaH0["deltahH2O"]
    h_co2 = h0_co2 + deltaH0["deltahCO2"]
    h_o2 = h0_o2 + deltaH0["deltahO2"]
    h_co = h0_co + deltaH0["deltahCO"]
    h_n2 = h0_n2 + deltaH0["deltahN2"]
    h_no = h0_no + deltaH0["deltahNO"]
    h_no2 = h0_no2 + deltaH0["deltahNO2"]
    h_oh = h0_oh + deltaH0["deltahOH"]
    h_h2 = h0_h2 + deltaH0["deltahH2"]
    h_h = h0_h + deltaH0["deltahH"]
    h_n = h0_n + deltaH0["deltahN"]
    h_o = h0_o + deltaH0["deltahO"]
    if combustivel=="H2":
        h_comb=h_h2
        h_h2=0
    else:
        h_comb = float(h0_comb + sp.integrate(Cp, (T, 298.15, i)))
    
    return h_h2o,h_co2,h_o2,h_co,h_n2,h_comb,h_no,h_no2,h_oh,h_h2,h_h,h_n,h_o

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
                    deltasComb=float(s0_comb + sp.integrate(Cp/T,(T,298.15,i)))
                
    return {"deltasComb":deltasComb,"deltasH2O":deltasH2O,"deltasCO2":deltasCO2,"deltasO2":deltasO2,"deltasCO":deltasCO,"deltasN2":deltasN2,"deltasNO":deltasNO,"deltasNO2":deltasNO2,"deltasOH":deltasOH,"deltasH2":deltasH2,"deltasH":deltasH,"deltasN":deltasN,"deltasO":deltasO}

def temperaturaAdiabatica(composicao,T0,T1):
    
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
        h_h2o,h_co2,h_o2,h_co,h_n2,h_comb,h_no,h_no2,h_oh,h_h2,h_h,h_n,h_o=interpolacaoH(T0)
        entalpy_products0=Ncomb*h_comb+NO2*h_o2+NH2O*h_h2o+NCO2*h_co2+NCO*h_co+NO*h_o+NN*h_n+NNO2*h_no2+NNO*h_no+NOH*h_oh+NH2*h_h2+NH*h_h+NN2*h_n2
        zero0=entalpy_products0-entalpy_reactants
        
        h_h2o,h_co2,h_o2,h_co,h_n2,h_comb,h_no,h_no2,h_oh,h_h2,h_h,h_n,h_o=interpolacaoH(T1)
        entalpy_products1=Ncomb*h_comb+NO2*h_o2+NH2O*h_h2o+NCO2*h_co2+NCO*h_co+NO*h_o+NN*h_n+NNO2*h_no2+NNO*h_no+NOH*h_oh+NH2*h_h2+NH*h_h+NN2*h_n2
        zero1=entalpy_products1-entalpy_reactants
        print(zero0,zero1)
        nextT=T1-zero1*(T1-T0)/(zero1-zero0)
        
        if abs(zero1)<0.5:
            break
        else:
            T0=T1
            T1=nextT
    return nextT

entalpy_reactants=(interpolacaoH(T0)[5]+OC*interpolacaoH(T0)[2]+3.76*OC*interpolacaoH(T0)[4])*mols
ent_prod=[entalpy_reactants]

bnds=((0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None))
cons=[{'type': 'eq', 'fun': restricaoC},{'type': 'eq', 'fun': restricaoO},{'type': 'eq', 'fun':restricaoH},{'type': 'eq','fun':restricaoN}]

estimative0 = (sys.float_info.min, sys.float_info.min, complete_combustion_H2O*mols, complete_combustion_CO2*mols, sys.float_info.min, sys.float_info.min,sys.float_info.min,sys.float_info.min,sys.float_info.min,sys.float_info.min,sys.float_info.min,sys.float_info.min,mols*3.76*OC)
#estimative0 = (sys.float_info.min, sys.float_info.min, complete_combustion_H2O*mols, 0, 0, sys.float_info.min, sys.float_info.min,sys.float_info.min,sys.float_info.min,sys.float_info.min,0,sys.float_info.min,mols*3.76*OC)
def composicaoAdiabatica(estimative,T0,P0,V):
    while True:
        
        P=P0*9.86923*10**-6
        
        sol = minimize(gibbs, estimative,args=(T0,P), method='SLSQP', bounds=bnds, constraints=cons)
        total_mols = sum(i for i in sol.x)
        
        nextTemp=float(temperaturaAdiabatica(estimative,T0,T0+1))
        print(nextTemp)
        P0=total_mols*Ru*nextTemp/V
        
        if abs(nextTemp-T0)<0.5:
            break
        else:
            T0=nextTemp
            estimative=sol.x
            
    return nextTemp,P0,sol.x,total_mols
    

for k in range(len(tk)):
   
    v=Vo.subs(teta,tk[k])
    Ta,Pa,moles,total_mols=composicaoAdiabatica(estimative0,Temp[k],Press[k],v)
    
    estimative0=moles
    
    Temp.append(Ta)
    Press.append(Pa)
    
    total_mass = sum([moles[x] * m[x] for x in range(13)])
        
    composicao_comb.append(moles[0] * m[0] / total_mass)
    composicao_o2.append(moles[1] * m[1] / total_mass)
    composicao_h2o.append(moles[2] * m[2] / total_mass)
    composicao_co2.append(moles[3] * m[3] / total_mass)
    composicao_co.append(moles[4] * m[4] / total_mass)
    composicao_o.append(moles[5] * m[5] / total_mass)
    composicao_n.append(moles[6] * m[6] / total_mass)
    composicao_no2.append(moles[7] * m[7] / total_mass)
    composicao_no.append(moles[8] * m[8] / total_mass)
    composicao_oh.append(moles[9] * m[9] / total_mass)
    composicao_h2.append(moles[10] * m[10] / total_mass)
    composicao_h.append(moles[11]* m[11] / total_mass)
    composicao_n2.append(moles[12] * m[12] / total_mass)
    
    h_h2o,h_co2,h_o2,h_co,h_n2,h_comb,h_no,h_no2,h_oh,h_h2,h_h,h_n,h_o=interpolacaoH(Ta) 
    
    entalpy_products=(moles[0]*h_comb+moles[1]*h_o2+moles[2]*h_h2o+moles[3]*h_co2+moles[4]*h_co+moles[5]*h_o+moles[6]*h_n+moles[7]*h_no2+moles[8]*h_no+moles[9]*h_oh+moles[10]*h_h2+moles[11]*h_h+moles[12]*h_n2)
    ent_prod.append(entalpy_products)
    
    print({'MASS_FRACTION':{'Combustivel': composicao_comb[k+1], 'O2':composicao_o2[k+1],
       'H2O':composicao_h2o[k+1], 'CO2': composicao_co2[k+1], 'CO': composicao_co[k+1],'O':composicao_o[k+1],'N':composicao_n[k+1],'NO2':composicao_no2[k+1],'NO':composicao_no[k+1],'OH':composicao_oh[k+1],'H2':composicao_h2[k+1],'H':composicao_h[k+1],'N2':composicao_n2[k+1]}})
    print({'MOLAR_FRACTION':{'Combustivel': moles[0] / total_mols, 'O2': moles[1] / total_mols, 'H2O': moles[2] / total_mols, 'CO2':moles[3]/total_mols, 'CO':moles[4]/total_mols,'O':moles[5]/ total_mols,'N':moles[6]/total_mols,'NO2':moles[7]/total_mols,'NO':moles[8]/ total_mols,'OH':moles[9]/ total_mols,'H2':moles[10]/ total_mols,'H':moles[11]/ total_mols,'N2':moles[12]/ total_mols}})

    print("TEMPERATURA ADIABATICA: "+str(Ta))
#plot das entalpias instantaneas
fig=plt.figure(figsize=(10,7))
ax=fig.add_subplot(1,1,1)

ax.scatter(Temp,list(map(lambda x:x*10**3,ent_prod)),color="y")
ax.legend(['Entalpias'])

fig2=plt.figure(figsize=(15,10))

plot=fig2.add_subplot(3,4,1)
plot.plot(Temp[1:] ,composicao_comb[1:],'o')
plot.legend(['Combustível'])

plot=fig2.add_subplot(3,4,2)
plot.plot(Temp[1:],composicao_o2[1:],'o')
plot.legend(['O2'])

plot=fig2.add_subplot(3,4,3)
plot.plot(Temp[1:],composicao_h2o[1:],'o')
plot.legend(['H2O'])

plot=fig2.add_subplot(3,4,4)
plot.plot(Temp[1:],composicao_co2[1:],'o')
plot.legend(['CO2'])

plot=fig2.add_subplot(3,4,5)
plot.plot(Temp[1:],composicao_co[1:],'o')
plot.legend(['CO'])

plot=fig2.add_subplot(3,4,6)
plot.plot(Temp[1:],composicao_o[1:],'o')
plot.legend(['O'])

plot=fig2.add_subplot(3,4,7)
plot.plot(Temp[1:],composicao_n[1:],'o')
plot.legend(['N'])

plot=fig2.add_subplot(3,4,8)
plot.plot(Temp[1:],composicao_no2[1:],'o')
plot.legend(['NO2'])

plot=fig2.add_subplot(3,4,9)
plot.plot(Temp[1:],composicao_no[1:],'o')
plot.legend(['NO'])

plot=fig2.add_subplot(3,4,10)
plot.plot(Temp[1:],composicao_oh[1:],'o')
plot.legend(['OH'])

plot=fig2.add_subplot(3,4,11)
plot.plot(Temp[1:],composicao_h2[1:],'o')
plot.legend(['H2'])

plot=fig2.add_subplot(3,4,12)
plot.plot(Temp[1:],composicao_h[1:],'o')
plot.legend(['H'])
#plt.subplot(3,4,13)
#plot=fig2.add_subplot(3,4,1)
#plot.plot(Temp,composicao_n2)
plt.show()
        