#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 23:20:30 2021

@author: luiza.maciel
"""


while True:
    combustivel=input("Com qual combustível deseja operar?  ")
    if combustivel=="gasolina":
        AC=9.6 #mols de oxigenio 
        ACm=12.99 #massica
        c=6.67
        hid=12.8
        ox=0.533
        concentracaoCO2=6.67
        concentracaoH2O=6.4
        concentracaoN2=3.76*AC
        PCI=39249*10**3 #kJ/kg
        a0=-8102.9675849653
        a1=7963.204888950
        a2=-2923.24238668569
        a3=509.144026930886
        a4=-42.2764373650608
        a5=1.35175803975684
        mol=((12*c+hid+16*ox)+AC*(32+3.76*28))/(1+AC+AC*3.76)  #massa molecular da mistura
        A=(5.7*10**11)
        E=(15098)
        m=0.25
        n=1.5
        break
    elif combustivel=="alcool":
        AC=3.2
        ACm=8.427
        c=2.15
        hid=6.62
        ox=1.23
        concentracaoCO2=2.15
        concentracaoH2O=3.31
        concentracaoN2=AC*3.76
        PCI=24804*10**3
        a0=-12482.8740213179
        a1=10263.4623453335
        a2=-3316.7850402598
        a3=526.309291795851
        a4=-40.8869367350809
        a5=1.24555084441151
        mol=((12*c+hid+16*ox)+AC*(32+3.76*28))/(1+AC+AC*3.76)
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
        concentracaoH2O=2.0722
        concentracaoCO2=1.0789
        concentracaoN2=AC*3.76+r
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