import sys

from scipy.optimize import minimize
import sympy as sp
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from math import log
from decimal import  Decimal
file2='./tabelaHS.xlsx'
file1='./dados.txt'
df2=pd.read_excel(file2)
arq=open(file1)
data=arq.readlines()
Temp=[]
#for t in data[3:303]:
#    Temp.append(float(t))
Temp=[300,600,900,1200,1500,2000]
carb=8*0.77+2*0.23
h=18*0.77+0.23*6
o=0.23
OC=9.6
h0_comb=0.77*(-208450)+0.23*(-235310)
concentracaoH2O=6.4
concentracaoCO2=6.67
h0_gas=-208450
h0_et=-235310
T=sp.symbols('T')
a0=-8102.9675849653
a1=7963.204888950
a2=-2923.24238668569
a3=509.144026930886
a4=-42.2764373650608
a5=1.3517580397568
P=1
mgas=114
met=46
mcomb=carb*12+o*16+h
m=[mcomb,32,18,44,28,16,14,46,30,17,2,1,28]
evolucao=[]
Cp=a0+a1*sp.log(T)+a2*sp.log(T)**2+a3*sp.log(T)**3+a4*sp.log(T)**4+a5*sp.log(T)**5  #J/molK
teta=T/1000
#Cp=a0+a1*sp.log(T)+a2*sp.log(T)**2+a3*sp.log(T)**3+a4*sp.log(T)**4+a5*sp.log(T)**5  #J/molK
#Cpgas=-0.056+6.75*teta-3.67*teta**2 +0.775*teta**3
Cpet=-0.2-4.65*teta-1.82*teta**2+0.03*teta**3
Ru=8.314
estimative=np.zeros(14)
def gibbs (composicoes):
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
    #Net=composicoes[5]
    Ntot=sum(composicoes)


    deltaG=Ncomb * (gt_comb + float((Decimal(Ncomb) * Decimal(P) / Decimal(Ntot)).ln())) + \
           NCO2 * (gt_co2 + float((Decimal(NCO2) * Decimal(P) / Decimal(Ntot)).ln())) + \
           NH2O * (gt_h2o + float((Decimal(NH2O) * Decimal(P) / Decimal(Ntot)).ln())) +\
           NO2 * (gt_o2 + float((Decimal(NO2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NO * (gt_o + float((Decimal(NO) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NN * (gt_n + float((Decimal(NN) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NNO2 * (gt_no2 + float((Decimal(NNO2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NNO * (gt_no + float((Decimal(NNO) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NOH * (gt_oh + float((Decimal(NOH) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NH2 * (gt_h2 + float((Decimal(NH2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NH * (gt_h + float((Decimal(NH) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NN2 * (gt_n2 + float((Decimal(NN2) * Decimal(P) / Decimal(Ntot)).ln()))+\
           NCO * (gt_co + float((Decimal(NCO) * Decimal(P) / Decimal(Ntot)).ln()))
    #print(deltaG)
    return (deltaG)


def restricaoC(composicoes):
    return carb-composicoes[0]*carb-composicoes[3]-composicoes[4]
def restricaoO(composicoes):
    return o+2*9.6-composicoes[0]*o-2*composicoes[1]-composicoes[2]-composicoes[4]-composicoes[5]-2*composicoes[7]-composicoes[8]-composicoes[9]-2*composicoes[3]
def restricaoH(composicoes):
    return h-h*composicoes[0]-2*composicoes[2]-composicoes[9]-2*composicoes[10]-composicoes[11]
def restricaoN(composicoes):
    return 9.6*3.76*2-composicoes[6]-composicoes[7]-composicoes[8]-2*composicoes[12]


#bnds=((sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None),(sys.float_info.min,None))
bnds=((0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None),(0,None))
cons=[{'type': 'eq', 'fun': restricaoC},{'type': 'eq', 'fun': restricaoO},{'type': 'eq', 'fun':restricaoH},{'type': 'eq','fun':restricaoN}]

for i in Temp:
    for k in range(36):
        if i>=df2['T'][k] and i<=df2['T'][k+1]:
            deltahH2O=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2O'][k]
            deltahCO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO2'][k]
            deltahO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO2'][k]
            deltahCO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hCO'][k]
            deltahNO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hNO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hNO'][k]
            deltahNO2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hNO2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hNO2'][k]
            deltahOH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hOH'][k]
            deltahH2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH2'][k]
            deltahH=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hH'][k]
            deltahN=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hN'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hN'][k]
            deltahN2=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hN2'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hN2'][k]
            deltahO=((i-df2['T'][k])/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k+1]+((df2['T'][k+1]-i)/(df2['T'][k+1]-df2['T'][k]))*df2['hO'][k]

        else:
            continue
    g0_h2o=-228590/(Ru*298.15)
    g0_co2=-394360/(Ru*298.15)
    g0_o2=0
    g0_co=-137150/(Ru*298.15)
    g0_no = 87600/(Ru*298.15)
    g0_no2 =51300/(Ru*298.15)
    g0_oh= 34280/(Ru*298.15)
    g0_h2 = 0
    g0_h = 203290/(Ru*298.15)
    g0_n = 455510/(Ru*298.15)
    g0_n2 = 0
    g0_o = 231770/(Ru*298.15)
    h0_h2o=-241820
    h0_o2=0
    h0_co2=-393520
    h0_co=-110530
    h0_no=90291
    h0_no2=33100
    h0_oh=39460
    h0_h2=0
    h0_h=218000
    h0_n=472650
    h0_n2=0
    h0_o=249190
    g0_gas=16530/(Ru*298.15)
    g0_et=(-168570)/(Ru*298.15)
    g0_comb=0.23*g0_et+0.77*g0_gas
    h_h2o=h0_h2o+deltahH2O
    h_co2=h0_co2+deltahCO2
    h_o2=h0_o2+deltahO2
    h_co=h0_co+deltahCO
    h_no = h0_no + deltahNO
    h_no2 = h0_no2 + deltahNO2
    h_oh = h0_oh + deltahOH
    h_h2 = h0_h2 + deltahH2
    h_h = h0_h + deltahH
    h_n = h0_n + deltahN
    h_n2 = h0_n2 + deltahN2
    h_o = h0_o + deltahO
    h_comb = h0_comb + sp.integrate(Cp, (T, 298.15, i))
   # h_gas=h0_gas+sp.integrate(Cpgas,(T,298.15,i))
    #h_et = h0_et + sp.integrate(Cpet, (T, 298.15, i))
    gt_h2o = g0_h2o - float(sp.integrate(h_h2o / (Ru * T ** 2), (T, 298.15, i)))
    gt_comb = g0_comb- float(sp.integrate(h_comb / (Ru * T ** 2), (T, 298.15, i)))
    #gt_gas = g0_gas - float(sp.integrate(h_gas / (Ru * T ** 2), (T, 298.15, i)))
    #gt_et = g0_et - float(sp.integrate(h_et / (Ru * T ** 2), (T, 298.15, i)))
    gt_co2 = g0_co2 - float(sp.integrate(h_co2 / (Ru * T ** 2), (T, 298.15, i)))
    gt_o2 = g0_o2 - float(sp.integrate(h_o2 / (Ru * T ** 2), (T, 298.15, i)))
    gt_co = g0_co - float(sp.integrate(h_co / (Ru * T ** 2), (T, 298.15, i)))
    gt_o = g0_o - float(sp.integrate(h_o / (Ru * T ** 2), (T, 298.15, i)))
    gt_n = g0_n - float(sp.integrate(h_n / (Ru * T ** 2), (T, 298.15, i)))
    gt_no2 = g0_no2 - float(sp.integrate(h_no2 / (Ru * T ** 2), (T, 298.15, i)))
    gt_no = g0_no - float(sp.integrate(h_no / (Ru * T ** 2), (T, 298.15, i)))
    gt_oh = g0_oh - float(sp.integrate(h_oh / (Ru * T ** 2), (T, 298.15, i)))
    gt_h2 = g0_h2 - float(sp.integrate(h_h2 / (Ru * T ** 2), (T, 298.15, i)))
    gt_h = g0_h - float(sp.integrate(h_h / (Ru * T ** 2), (T, 298.15, i)))
    gt_n2 = g0_n2 - float(sp.integrate(h_n2 / (Ru * T ** 2), (T, 298.15, i)))
    #est=(0,0,concentracaoH2O,concentracaoCO2,0,0,0,0,0,0,0,0,9.6*3.76)
    est = (sys.float_info.min, sys.float_info.min, concentracaoH2O, concentracaoCO2, sys.float_info.min, sys.float_info.min, sys.float_info.min, sys.float_info.min, sys.float_info.min, sys.float_info.min, sys.float_info.min, sys.float_info.min, 9.6 * 3.76)
    #est=(1,1,1,1,1,1,1,1,1,1,1,1,1)
    sol=minimize(gibbs,est ,method='SLSQP',bounds=bnds,constraints=cons)

    total_mass=sum([sol.x[x]*m[x]for x in range(13)])
    tot=sum(i for i in sol.x)
    evolucao.append(sol.x[0]/tot)
    print({'Combustivel':sol.x[0]*m[0]/total_mass,'O2':sol.x[1]*m[1]/total_mass,'H2O':sol.x[2]*m[2]/total_mass,'CO2': sol.x[3]*m[3]/total_mass, 'CO':sol.x[4]*m[4]/total_mass,'O':sol.x[5]*m[5]/total_mass,'H2':sol.x[11]*m[11]/total_mass,'NO2':sol.x[8]*m[8]/total_mass})
    print({'combsutivel':sol.x[0]/tot,'O2':sol.x[1]/tot,'H2O':sol.x[2]/tot,'CO2': sol.x[3]/tot, 'CO':sol.x[4]/tot,'O':sol.x[5]/tot,'H2':sol.x[11]/tot,'NO2':sol.x[7]/tot})

plt.plot(Temp,evolucao)

plt.show()




