# Importing required modules
import streamlit as st
import sympy as sym
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Constants (unchanged)
L = 2.99e8               # Speed of Light in m/s
r = 2.817e-15            # Classical Radius of Electron
R = 10                   # Lorentz Function Gamma 1
Q = 10**5                # Lorentz Function Gamma 2
m = 9e-31                # Mass of Electron
s = 2.3676e12            # Constant C3
h = 6.626e-34            # Planck's Constant
k = 1.38e-23             # Boltzmann's Constant

# Sidebar Inputs (unchanged)
st.sidebar.write("Input values")
p = st.sidebar.number_input("Enter value for Electron Index: ", value=2.5)
B = st.sidebar.number_input("Enter value for Magnetic Field in T: ", value=2.72)
q = st.sidebar.number_input("Enter value for L: ", value=2.72)
alpha = st.sidebar.number_input("Enter value for Seed Photon Index:", value=2.72)
F = st.sidebar.number_input("Enter value for Flux Density in W.m^-2.Hz^-1:", value=2.72)
v = st.sidebar.number_input("Enter value for Reference Frequency in Hz:", value=2.72)
T = st.sidebar.number_input("Enter value for Temperature:", value=None)
lower_E = st.sidebar.number_input("Enter value for lower limit of epsilon: ", value=1)
upper_E = st.sidebar.number_input("Enter value for upper limit of epsilon: ", value=100)

# User Defined Functions (updating with Trapezoidal Rule integration)

# Trapezoidal Rule Function
def trapezoidal_rule(E, V, p):
    """
    Applies Trapezoidal Rule to the integral of V(E)*E^((p-1)/2) over E.
    Assumes E is equally spaced.
    """
    n = len(E) - 1
    h = E[1] - E[0]
    
    # Check if epsilon values are equally spaced
    if not np.allclose(np.diff(E), h, rtol=1e-5, atol=1e-8):
        raise ValueError("Epsilon values must be equally spaced.")
    
    # Transform the function as V(E) * E^((p-1)/2)
    transformed = (V * E**((p - 1) / 2)).to_numpy()
    
    # Apply Trapezoidal Rule
    result = transformed[0] + transformed[-1] + 2 * np.sum(transformed[1:n])
    
    # Compute the integral approximation
    P2 = (h / 2) * result
    return P2

# Function to calculate int2 for the logarithmic difference (unchanged)
def int2(E, p, alpha):
    Ei = E.iloc[0]
    Ef = E.iloc[-1]
   
    if alpha == (p-1)/2:
        if Ei > 0:
            diff = np.log(Ef) - np.log(Ei)
        else:
            raise ValueError("Logarithm undefined for non-positive values.")
    else:
        diff = ((Ef**(((p-1)/2)-alpha))/(((p-1)/2)-alpha)) - ((Ei**(((p-1)/2)-alpha))/(((p-1)/2)-alpha))
        
    return diff

# Data handling and integration updates
def createdata(file):
    df = pd.read_csv(file, delim_whitespace=True)
    if 'Epsilon' not in df.columns or 'V_Epsilon' not in df.columns:
        st.error("Required columns 'Epsilon' or 'V_Epsilon' are missing in the dataset.")
        return
    E = df['Epsilon']
    V = df['V_Epsilon']
    return E, V

# File uploader logic (unchanged)
st.sidebar.write("Upload a text file containing data of curve of V(Epsilon) against Epsilon")
uploaded_file = st.sidebar.file_uploader("Choose a text file", type="txt")
st.sidebar.write("**columns should be named 'Epsilon', 'V_Epsilon'")

if uploaded_file is None:
    st.sidebar.write("Using sample dataset. Upload file and provide values to process other datasets.")
    E, V = createdata('Sample.txt')
   
if uploaded_file is not None:
    E, V = createdata(uploaded_file)

# Calculate the integrals using Trapezoidal Rule
P2_trap = trapezoidal_rule(E, V, p)
diff = int2(E, p, alpha)

# Calculate constants and final values for Volume Emissivity (unchanged)
o = (F * (v**alpha)) / (L * (h**(1-alpha)))
P1 = Const1(p)
P1T = Const1T(p, T) if T is not None else None

# Final result calculation
if T is None:
    Constt1 = P1 * P2_trap
else:
    Constt1 = P1T * P2_trap

Constt2 = P1 * o * diff

# Creating lists for plotting (unchanged)
final = [Constt1 * (i**(-((p-1)/2))) * 1.6e-16 for i in liste]
final2 = [Constt2 * (i**(-((p-1)/2))) * 1.6e-16 for i in liste]

# Data for display (updated with new results)
data = pd.DataFrame({'Epsilon': liste, 'Volume Emissivity': final})
data2 = pd.DataFrame({'Epsilon': liste, 'Volume Emissivity': final2})

data2["Epsilon"] = data["Epsilon"].apply(lambda x: '{:.6e}'.format(x))
data2["Volume Emissivity"] = data["Volume Emissivity"].apply(lambda x: '{:.6e}'.format(x))
data["Epsilon"] = data["Epsilon"].apply(lambda x: '{:.6e}'.format(x))
data["Volume Emissivity"] = data["Volume Emissivity"].apply(lambda x: '{:.6e}'.format(x))

# Plotting functions (unchanged)
def plot_it(liste, final, x_label, y_label, title):
    plt.figure(figsize=(10, 6))
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.plot(liste, final, color='grey', alpha=0.5)
    plt.title(title)
    logscale = True #st.toggle("Show Graph in logscale", value=True)
    if logscale:
        plt.xscale('log')
        plt.yscale('log')
    plt.grid(True)
    plt.legend()
    st.pyplot(plt)

# Streamlit Layout (unchanged)
st.title("Inverse Compton Spectra for Single Scattering")

st.header("Power Law Distribution")
st.write("")
st.dataframe(data2, use_container_width=True)
plot_it(liste, final2, 'Epsilon1', 'Volume Emissivity', 'Inverse Compton Result')
st.latex(r' \frac{dE}{dVdtd\epsilon_{1}} = \pi cr_{0}^2 CA(p)\epsilon_{1}^\frac{{-\left(p-1 \right)}}{2} \int_{}^{}d\epsilon \epsilon^\frac{{\left(p-1 \right)}}{2} v(\epsilon) ')
'Where'
st.latex(r' A(p) = 2^{p+3} \frac{p^2 + 4p + 11}{(p+3)^2 (p+5)(p+1)} ')
'And'
st.latex(r' C = N_{0}(m_{e}c^2)^{-p} ')

st.write("")
st.header("Black Body Condition")
st.write("")
st.dataframe(data, use_container_width=True)
plot_it(liste, final, 'Epsilon1', 'Volume Emissivity', 'Inverse Compton Result')
'Assuming value of T is not zero'
st.latex(r' \frac{dE}{dVdtd\epsilon_{1}} = \frac{C8\pi^2r_{0}^2}{h^3c^2} \left(kT \right)^\frac{\left( p+


st.markdown("""
<div style='text-align: right;'>
    <p><strong>By Garv Trivedi</strong></p>
    <p><strong>under guidance of Dr. C. Konar</strong></p>
</div>
""", unsafe_allow_html=True)

