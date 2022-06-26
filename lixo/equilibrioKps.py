# -*- coding: utf-8 -*-
"""
Created on Wed Apr 14 22:04:40 2021

@author: lumac
"""
import sys

import pandas as pd
import sympy as sp
from scipy.optimize import *
from decimal import *
file='./tabelaKps.xlsx'
file2='./tabelaHS.xlsx'
df=pd.read_excel(file,dtype={'Kp1':float,'Kp2':float,'Kp3':float,'Kp4':float,'Kp5':float,'Kp6':float,'Kp7':float,'Kp8':float})
df2=pd.read_excel(file2)


Temp=[300,500,450,1100]
Temp=input('temperaturas: ').split(',')

P=[1,3,5,4]
P=Decimal(1)
c=Decimal(6.67)
h=Decimal(12.8)
o=Decimal(0.533)
OC=Decimal(9.6)
h0=0.77*(-208600)+0.23*(-235000)
s0=0.77*(466.514)+0.23*(282.444)
T=sp.symbols('T')
a0=-8102.9675849653
a1=7963.204888950
a2=-2923.24238668569
a3=509.144026930886
a4=-42.2764373650608
a5=1.35175803975684
Cp=a0+a1*sp.log(T)+a2*sp.log(T)**2+a3*sp.log(T)**3+a4*sp.log(T)**4+a5*sp.log(T)**5
gcomb=h0+sp.integrate(Cp,(T,298.15,Temp[0]))-T*(s0+sp.integrate(Cp/T,(T,298.15,Temp[0])))
Ru=8.314
concentracaoH2O=6.4
concentracaoCO2=6.67
concentracaoO2=9.60
#calcula 
def composicoes(lista,P):
    Ncomb=Decimal(lista[0])
    NH2O=Decimal(lista[1])
    NH2=Decimal(lista[2])
    NO2=Decimal(lista[3])
    NO=Decimal(lista[4])
    NOH=Decimal(lista[5])
    NN2=Decimal(lista[6])
    NN=Decimal(lista[7])
    NH=Decimal(lista[8])
    NCO=Decimal(lista[9])
    NCO2=Decimal(lista[10])
    NNO=Decimal(lista[11])
    NNO2=Decimal(lista[12])
    Ntot=NH2O+NH2+NO2+NO+NOH+NN2+NH+NCO+NCO2+NNO+NNO2+NN+Ncomb
    meio=Decimal(0.5)
    eq1=((NH**2/NH2)*(P/Ntot))-Decimal(1e1)**(Kp1)
    eq2=((NO**2/NO2)*(P/Ntot))-Decimal(1e1)**(Kp2)
    eq3=((NN**2/NN2)*(P/Ntot))-Decimal(1e1)**(Kp3)
    eq4=((NH2*NO2**meio/NH2O)*(P/Ntot)**meio)-Decimal(1e1)**(Kp4)
    eq5=((NH2**meio*NOH/NH2O)*(P/Ntot)**meio)-Decimal(1e1)**(Kp5)
    eq6=((NCO*NO2**meio/NCO2)*(P/Ntot)**meio)-Decimal(1e1)**(Kp6)
    eq7=(NNO/(NN2**meio*NO2**meio))-Decimal(1e1)**(Kp7)
    eq8=((NNO2/(NN2**meio*NO2))*(P/Ntot)**-meio)-Decimal(1e1)**(Kp8)
    eq9=((NH2O**Decimal(concentracaoH2O)*NCO2**Decimal(concentracaoCO2)/(NO2**Decimal(concentracaoO2)*Ncomb))*(P/Ntot)**(Decimal(concentracaoH2O)+Decimal(concentracaoCO2)-Decimal(concentracaoO2)-1))-Decimal(1e1)**(Kp9)
    eq10=NCO+NCO2+Ncomb*c-c
    eq11=2*NH2O+2*NH2+NOH+NH+Ncomb*h-h
    eq12=2*NN2+NNO+NNO2-OC*Decimal(3.76)*2
    eq13=NH2O+2*NO2+NO+NOH+NCO+2*NCO2+NNO+2*NNO2+Ncomb*o-o-2*OC

    return eq1,eq2,eq3,eq4,eq5,eq6,eq7,eq8,eq9,eq10,eq11,eq12,eq13
estimative=input('digite uma estimativa inicial para as composicoes: ')
if estimative =='':
    #estimative=[1,1,1,1,1,1,1,1,0,1,1,1,1]
    estimative=[sys.float_info.min,concentracaoH2O,sys.float_info.min,sys.float_info.min,sys.float_info.min,sys.float_info.min,3.76*9.6,sys.float_info.min,sys.float_info.min,sys.float_info.min,concentracaoCO2,sys.float_info.min,sys.float_info.min]
else:
    estimative=estimative.split(',')
evolucao=[]

for t in Temp:
    i=float(t)
    for j in range(22):
        if i >= df['T'][j] and i<=df['T'][j+1]:
            Kp1=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp1'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp1'][j])
            Kp2=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp2'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp2'][j])
            Kp3=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp3'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp3'][j])
            Kp4=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp4'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp4'][j])
            Kp5=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp5'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp5'][j])
            Kp6=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp6'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp6'][j])
            Kp7=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp7'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp7'][j])
            Kp8=Decimal(((i-df['T'][j])/(df['T'][j+1]-df['T'][j]))*df['Kp8'][j+1]+((df['T'][j+1]-i)/(df['T'][j+1]-df['T'][j]))*df['Kp8'][j])
            
        else:
            continue
    for k in range(36):
        if i>=df2['T'][k] and i<=df2['T'][k+1]:
            deltahH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k]
            deltasH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sH2O'][k]
            deltahCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k]
            deltasCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sCO2'][k]
            deltahO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k]
            deltasO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['sO2'][k]

        else:
            continue

    h0_h2o = -241820
    h0_o2 = 0
    h0_co2 = -393520
    h_comb = h0 + float(sp.integrate(Cp, (T, 298.15, i)))
    s_comb=s0+float(sp.integrate(Cp/T,(T,298.15,i)))
    h_h2o = h0_h2o + deltahH2O
    h_co2 = h0_co2 + deltahCO2
    h_o2 = h0_o2 + deltahO2
    g_comb=h_comb-i*s_comb
    g_h2o=h_h2o-i*deltasH2O
    g_co2=h_co2-i*deltasCO2
    g_o2=h_o2-i*deltasO2
    delta_G=Decimal(g_co2*concentracaoCO2+g_h2o*concentracaoH2O-g_o2*concentracaoO2-g_comb)
    Kp9=(-delta_G/(Decimal(Ru)*Decimal(i))).ln()
    print(Kp9)
    x0=[]
    for e in estimative:
        x0.append(float(e))
    sol=fsolve(composicoes,x0,args=1)
    print (sol)#para a temperatura
    evolucao.append(sol)
    print (composicoes(sol,1))

# g0_h2o = -228590
# g0_co2 = -394360
# g0_o2 = 0
  # delta0_G = g0_co2 * concentracaoCO2 + g0_h2o * concentracaoH2O - g0_o2 * concentracaoO2 - g0_comb
    # deltaH0 = h_co2 * concentracaoCO2 + h_h2o * concentracaoH2O - h_o2 * concentracaoO2 - h_comb
    # deltaH = sp.integrate(deltaH0 / (Ru * T ** 2), (T, 298.15, i))  #vanthoff
    # delta_GRT = deltaH + delta0_G / (Ru * 298.15)
