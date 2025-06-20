import streamlit as st
import sympy as sym
st.header ("Wow")
#constants
h = 6.626 * (10**(-34))
L = 2.99 * (10**8)
k = 1.38
r = 2.817 * (10**(-15))

print("input values")
p = st.number_input("Enter value for Momentum: ",value=2.5)
E = st.number_input("Enter value for Energy of the Seed Photon: ",value=2.72)
T = st.number_input("Enter value for Temperature in K: ",value=2.72)
c = st.number_input("Enter value for Normalization of electron energy distribution: ",value=2.72)

def vol_emmissivity(p,E,T,c):
    a = p + 3
    #print(f"a={a}")                                                        
    b = (p + 5) / 2
    #print(f"b={b}")                                                        
    x = (p - 1) / 2
    #print(f"x={x}")                                                        
    A = ((p**2) + 4*p + 11) / ((a**2)*(2*b)*(p + 1))
    #print(f"A={A}")
    g = sym.gamma(b)
    #print(f"g={g}")
    Z = sym.zeta(b)
    #print(f"Z={Z}")
    F = A*g*Z
    #print(f"F={F}")
    i = 3.141592653589
    d = (8*(i**2)*(r**2)) / ((h**3)*(L**2))
    #print(f"d={d}")
    o = (k*T)**b
    #print(f"o={o:e}")

    X = (c*d*o*F*(E*(-x)))*((1e-23)**b)
    exp = X *(1.6*(10**(-16)))
    return exp

value=vol_emmissivity(p,E,T,c)
st.info(f"Volume Emmissivity is {value:e} Js-1KeV-1K-1")
