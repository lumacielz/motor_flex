#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Feb  5 11:34:38 2022

@author: luiza.maciel
"""


from math import radians,degrees,sin,cos,sqrt,log,exp
import numpy as np
import sympy as sp
from scipy import signal
from scipy.optimize import curve_fit


#dados da geometria do motor
dados_motor=input('Dados do Motor: ').split(',')
#motor - tadeu
#dados_motor=11,0.144,0.081,0.0864,4

r=float(dados_motor[0])
L=float(dados_motor[1])
D=float(dados_motor[2])
S=float(dados_motor[3])
cilindros=float(dados_motor[4])

avanco=float(input("Ângulo de avanço: "))
if avanco<=45 and avanco>36:
    fcor=(avanco-36)*((1.77-1.67)/(45-36))+1.67
elif avanco<=36 and avanco >27:
    fcor=(avanco-27)*((1.67-1.68)/(36-27))+1.68
elif avanco<=27 and avanco >=18:
    fcor=(avanco-18)*((1.68-1.45)/(27-18))+1.45
else:
    fcor=1

print('fcor = '+str(fcor))
teta0d=float(input("ângulo de fechamento da válvula de admissão [em graus]: "))
teta_combd=float(input("ângulo de início da combustão[em graus]: "))
delta_tetad=float(input("duração da combustão[em graus]: "))
tetafd=float(input("ângulo de abertura da válvula de escape[em graus] "))

teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)

#condicoes de operacao
Var=float(input("Vazão de ar[kg/h]: "))
P0=float(input("Pressão no ângulo de fechamento: "))*10**3
Tp=float(input("Temperatura das paredes: "))
rotacao=float(input("rotação do motor: "))

R=S/2 #raio do virabrequim
N=rotacao*2*np.pi/60 #rad/s
Vd=np.pi*D**2*R/2 #volume deslocado m^3
vp=2*S*N #velocidade do pistão m/s

Ru=8.314 #kJ/kmolK
M_air = 28.962


def stringToFloat(inputLista):
    lista=[]
    for i in inputLista.split(","):
        lista.append(float(i))
    return lista

def composicaoCombustivel(c):
    coef=c.split("C")[1].split("H")
    coef[1:]=(coef[1].split("O"))
    return list(map(lambda x: float(x),coef))

def Cps(Ts):
    CpH2O = 0.14904476*log(Ts)**5 -5.60452380*log(Ts)**4 + 81.54220310*log(Ts)**3 -573.43988184*log(Ts)**2 + 1954.64102748*log(Ts) -2558.72544088
    CpC2H5OH = 1.40634*log(Ts)**5 -45.60008*log(Ts)**4 + 578.32500*log(Ts)**3 -3579.01508*log(Ts)**2+ 10834.12151 *log(Ts)-12834.47773
    CpC8H18 = 4.9695987*log(Ts)**5 -157.61132*log(Ts)**4 +1966.430702 *log(Ts)**3 -12036.09466*log(Ts)**2 +36241.19904*log(Ts)-43029.69896
    return CpC8H18*oc+et*CpC2H5OH+h2o*CpH2O

def fitLog(Ts, a0,a1,a2,a3,a4,a5):
    return a0+a1*np.log(Ts)+a2*np.log(Ts)**2+a3*np.log(Ts)**3+a4*np.log(Ts)**4+a5*np.log(Ts)**5

tempRange = np.linspace(298.15,3000,1000)

while True:
    combustivel=input("Com qual combustível deseja operar?  ")
    if combustivel=="gasolina":
        ocV = float(input("%V octano: "))
        etV = float(input("%V etanol: "))
        
        n_oc = ocV*703/(114.232)
        n_et = etV*789/(46.069)
        oc = n_oc/(n_oc+n_et)
        et = n_et/(n_oc+n_et)
 
        h2o=0
        c,h,o,n=8*oc+2*et,18*oc+et*6,et,0
        
        cps = [Cps(t) for t in tempRange]
        a0,a1,a2,a3,a4,a5=curve_fit(fitLog,tempRange,cps)[0]

        h0_comb=oc*(-208450)+et*(-235310)
        s0_comb=oc*(466.514)+et*(282.444)

        PCI = oc*44430*10**3 + et*26810*10**3
        
        Ae=(5.7*10**11) #octano
        Ea=15098
        me=0.25
        ne=1.5
        
        break
        
    elif combustivel=="alcool":
#        et=1
#        h2o=0
        et = 0.938756707502278
        h2o = 0.06124329249772197
        oc=0
        c,h,o,n = et*2,6*et+2*h2o,h2o+et,0
        
        cps = [Cps(t) for t in tempRange]
        a0,a1,a2,a3,a4,a5=curve_fit(fitLog,tempRange,cps)[0]
    
        h0_comb=et*(-235310)+h2o*(-241826)
        s0_comb=et*(282.444)+h2o*(188.835)

        PCI= 24804*10**3
        
        Ae=(1.8*10**12) #etanol
        Ea=15098
        me=0.15
        ne=1.6
        
        break
    
    elif combustivel=="GNV":
        
        met=0.92285
        et=5.455*10**-2
        prop=1.115*10**-2
        but=0.125*10**-2
        co2=0.255*10**-2
        n=0.765*10**-2
       # met,et,prop,but,co2,n=1,0,0,0,0,0
        
        c=met+2*et+3*prop+4*but+co2
        h=met*4+et*6+prop*8+but*10
        o=co2*2
    
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
        h0_comb=met*h0_met+et*h0_et+prop*h0_prop+but*h0_but+co2*h0_co2+n*h0_n2
        s0_comb=met*s0_met+et*s0_et+prop*s0_prop+but*s0_but+co2*s0_co2+n*s0_n2
      
        PCI = 48737*10**3
        
        Ae=1.59*10**13
        beta=0
        Ea=47.8/(1.987*10**-3)
        me=0.7
        ne=0.8
        
        break 
    
    elif combustivel == "H2":
        c,h,o,n=0,2,0,0
        mcomb=2
        a5 = -0.1571988
        a4 = 5.07583945
        a3 = -64.62974617
        a2 = 406.63764375
        a1 = -1265.40995607
        a0 = 1586.40614067
        
        PCI = 120000*10**3
        break
    elif combustivel == "outro":
        molecularComposition = input("digite a composição do combustível na forma CnHnOn: ")
        c,h,o=composicaoCombustivel(molecularComposition),0
        n=0
        mcomb=c*12+h+o*16
        
        PCI = float(input("Poder Calorífico Inferior: "))
        
        h0_comb=float(input("entalpia de formação do combustivel: "))
        s0_comb=float(input("entropia de formação do combustivel: "))
        
        a0,a1,a2,a3,a4,a5=stringToFloat(input("coeficientes da equação logaritimica de Cp em funcao da temperatura: "))
     
    else:
        print("digite um combustível válido!!!")

#massa molecular do combustivel        
M_c=(12.0107*c+h*1.008+16*o+28*n)

#quantidades estequimetricas de produtos
stoichometricH2O = h/2
stoichometricCO2 = c
OC = (2*stoichometricCO2 + stoichometricH2O - o)/2 #mols de oxigenio
stoichometricN2 = 3.76*OC+n/2

#TODO verificar massa do ar
AC = OC*(32 + 3.76 * 28.16) / M_c  #razao ar combustivel
M_m=(M_c+OC*(32+3.76*28.16))/(1+OC+OC*3.76)  #massa molecular da mistura

#massa de ar admitida
m_air=Var/(0.5*3600*cilindros*(N/(2*np.pi))) #rot por s
#massa de mistura admitida
m_m=m_air/AC+m_air 
#massa de combustivel admitida
m_c=m_air/AC
mols_c = m_c/M_c

#fracao massica inicial
Y0 = m_c/m_m
Rg=Ru/M_m #constante dos gases

#DADOS PARA WIEBE!
coefwiebe=input('Coeficientes da equação de Wiebe m , a: ').split(',')
m=float(coefwiebe[0])
a=float(coefwiebe[1])
#m,a=2,5
ef=0.87
Qtot=ef*m_c*PCI #kJ

#DADOS PARA EQUILIBRIO
#mols de oxigenio e nitrogenio admitidos
O=OC*mols_c
NI=3.76*O
X0 = mols_c/(mols_c+O+NI)

#composicao molar inicial
reactants_composition = np.asarray([mols_c,O,0,0,0,0,0,0,0,0,0,0,NI])
#comb,o2,h2o,co2,co,o,n,no2,no,oh,h2,h,n2
estimative = np.asarray([10**-15,10**-15,stoichometricH2O*mols_c,stoichometricCO2*mols_c,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,NI])
#estimative = np.asarray([mols_c,OC*mols_c,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,10**-15,(NI+n*mols_c)])
Molecular_Mass=[M_c,31.999,18.015,44.01,28.01,16,14.007,46.005,30.006,17.007,2.016,1.008,28.013]

#TODO verificar se e necessario
Ts = sp.symbols('Ts')
Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5

#massa para cinetica
Molecular_Mass_Kinects=[M_c,32,28,44,18,28]
comp=[reactants_composition[0]]
