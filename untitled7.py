#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 08:04:17 2022

@author: luiza.maciel
"""


import sympy as sp

a0=-8102.9675849653
a1=7963.204888950
a2=-2923.24238668569
a3=509.144026930886
a4=-42.2764373650608
a5=1.35175803975684

Ts = sp.symbols('Ts')
h0_comb=0.5462*(-208450)+0.4538*(-235310)
Cp=a0+a1*sp.log(Ts)+a2*sp.log(Ts)**2+a3*sp.log(Ts)**3+a4*sp.log(Ts)**4+a5*sp.log(Ts)**5
deltahc = sp.integrate(Cp,(Ts,298.15,700))
hc=deltahc-208600
ho2=12499
hh2o=14190-241826
hco2=17754-393522
hco=110527+12021
ho=8570+249170
hn=8353+472680
hno2=33100+17250
hno=12308+90291
hoh=11902+38987
hh2=11730
hh=217999+8353
hn2=11937

hs=[float(hc),ho2,hh2o,hco2,hco,ho,hn,hno2,hno,hoh,hh2,hh,hn2]

c0=[2.647873991049049e-07, 2.168317532530156e-06, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 8.181062050236279e-06]
c1=[3.01371229e-13, 7.62521803e+00, 4.32054021e-01, 3.63224922e-01,1.45828882e-12, 5.32780602e-13, 9.27388920e-13, 4.01054587e-05,1.45028324e-05, 6.06456223e-13, 1.58213009e-12, 1.46997020e-13,3.07902367e+01]
delta=0
for c,h in zip(c1,hs):
    delta += c*h
print(delta)
