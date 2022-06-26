#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Nov 21 13:52:01 2021

@author: luiza.maciel
"""


from math import radians,exp
import matplotlib.pyplot as plt
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp,odeint


print("\n------------------------------ANALISE DE MOTOR FLEX------------------------------\n")

Ru=8.314
#dados_motor=input('Dados do Motor: ').split(',')
dados_motor=[11,0.144,0.081,0.086,4]
r=float(dados_motor[0])
L=float(dados_motor[1])
D=float(dados_motor[2])
S=float(dados_motor[3])
cilindros=float(dados_motor[4])

#avanco=float(input("Ângulo de avanço: "))
avanco=21
if avanco<=45 and avanco>36:
    fcor=(avanco-36)*((1.77-1.67)/(45-36))+1.67
elif avanco<=36 and avanco >27:
    fcor=(avanco-27)*((1.67-1.68)/(36-27))+1.68
elif avanco<=27 and avanco >=18:
    fcor=(avanco-18)*((1.68-1.45)/(27-18))+1.45
print('fcor= '+ str(fcor))

#teta0d=float(input("ângulo de fechamento da válvula de admissão [em graus]: "))
#teta_combd=float(input("ângulo de início da combustão[em graus]: "))
#delta_tetad=float(input("duração da combustão[em graus]: "))
#tetafd=float(input("ângulo de abertura da válvula de escape[em graus] "))
teta0d= -164
tetafd=146
teta_combd=-4
delta_tetad=38
teta0=radians(teta0d)
tetaf=radians(tetafd)
teta_comb=radians(teta_combd)
delta_teta=radians(delta_tetad)

rotacao=2498
P0=67.49*10**3
Var=89.25
Tp=373
#Var=float(input("Vazão de ar[kg/h]: "))
#P0=float(input("Pressão no ângulo de fechamento: "))*10**3
#Tp=float(input("Temperatura das paredes: ")) #fazer input
 #kJ/kmolK
#rotacao=float(input("rotação do motor: "))
N=rotacao*2*3.14/60 #rad/s
R=S/2 #raio do virabrequim
Vd=3.14*D**2*R/2 #volume deslocado m^3
vp=2*S*N #velocidade do pistão m/s

while True:
    combustivel=input("Com qual combustível deseja operar?  ")
    if combustivel=="gasolina":
        AC=9.6 #mols de oxigenio 
        ACm=12.99 #massica
        c=6.67
        hid=12.8
        ox=0.533
        stoichometricCO2=6.67
        stoichometricH2O=6.4
        stoichometricN2=3.76*AC
        PCI=39249*10**3 #kJ/kg
        a0=-8102.9675849653
        a1=7963.204888950
        a2=-2923.24238668569
        a3=509.144026930886
        a4=-42.2764373650608
        a5=1.35175803975684
        mol=((12*c+hid+16*ox)+AC*(32+3.76*28))/(1+AC+AC*3.76)  #massa molecular da mistura
        Ae=(5.7*10**11) #octano
        Ea=15098
        me=0.25
        ne=1.5
        break
    elif combustivel=="alcool":
        AC=3.2
        ACm=8.427
        c=2.15
        hid=6.62
        ox=1.23
        stoichometricCO2=2.15
        stoichometricH2O=3.31
        stoichometricN2=AC*3.76
        PCI=24804*10**3
        a0=-12482.8740213179
        a1=10263.4623453335
        a2=-3316.7850402598
        a3=526.309291795851
        a4=-40.8869367350809
        a5=1.24555084441151
        mol=((12*c+hid+16*ox)+AC*(32+3.76*28))/(1+AC+AC*3.76) #massa molecular da mistura
        Ae=(1.8*10**12) #etanol
        Ea=15098
        me=0.15
        ne=1.6
        break
    elif combustivel=="GNV":
        AC=2.112 #molar
        ACm=16.64
        met=0.92285
        #92.285 em volume
        et=5.455*10**-2
        prop=1.115*10**-2
        but=0.125*10**-2
        co2=0.255*10**-2
        r=0.765*10**-2
        stoichometricH2O=2.0722
        stoichometricCO2=1.0789
        stoichometricN2=AC*3.76+r
        PCI=48737*10**3
        a0=-12713.9307245771
        a1=10272.7734235011
        a2=-3251.2253041433
        a3=504.476942937525
        a4=-38.3629908265519
        a5=1.14655193729902
        mol=((met*16+et*30+prop*44+but*58+co2*44+28*r)+AC*(32+3.76*28))/(1+AC+AC*3.76)
        break
    else:
        print("digite um combustível válido!!!")

m_ar=Var/(0.5*3600*cilindros*(N/(2*3.14))) #rot por s
m_m=m_ar/ACm+m_ar #massa da mistura admitida
m_c=(m_m/(1+ACm)) 
m_c=m_ar/ACm
print(m_m,m_c)  
m_cm=hid+c*12+ox*16
print(m_c/m_cm,m_m/mol)
m0=1/((1+AC)+3.76*AC)*m_cm    

Ts=sp.symbols('Ts')
Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5
CpO=10228.3426-7184.92333*sp.log(Ts)+2010.86808*sp.log(Ts)**2-279.69496*sp.log(Ts)**3+19.34823*sp.log(Ts)**4-0.53257*sp.log(Ts)**5
CpN=0.54497*sp.log(Ts)**5-18.69984*sp.log(Ts)**4+254.29554*sp.log(Ts)**3-1712.17390*sp.log(Ts)**2+5708.38047*sp.log(Ts)-7513.3642
CpCO2=0.20754*sp.log(Ts)**5-6.43522*sp.log(Ts)**4+77.54809*sp.log(Ts)**3-452.81197*sp.log(Ts)**2+1288.4677*sp.log(Ts)-1412.36785
CpH2O=0.64541*sp.log(Ts)**5-23.54277*sp.log(Ts)**4+339.33662*sp.log(Ts)**3-2414.77575*sp.log(Ts)**2+8490.5218*sp.log(Ts)-11780.76495

KReag=((Cp/(Cp-Ru))+AC*(CpO/(CpO-Ru))+3.76*AC*(CpN/(CpN-Ru)))/(1+AC+3.76*AC)
KProd=((CpCO2/(CpCO2-Ru))*stoichometricCO2+(CpH2O/(CpH2O-Ru))*stoichometricH2O+(CpN/(CpN-Ru))*AC*3.76)/(stoichometricCO2+stoichometricH2O+stoichometricN2)

ef=0.87
Qtot=ef*m_c*PCI #kJ
Rg=Ru/mol #constante dos gases
