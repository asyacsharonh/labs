#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np

data = np.loadtxt('INSERT FILE NAME', delimiter='')

#VALUES, denoted as (variable) and d(variable) for the uncertainty
V0_exp = data[:,0] #from csv data collected in lab
dV0 = 

d = #from measured thickness of plastic separation  
dd = 

n = 18.5e-6  #NOTE: might change due to temp, this value is for 20deg
dn =

v0_exp = data[:,1] #from csv data
dv0 =

p = #density of oil drop
dp =

sigma = #density of air
dsigma =

g = 9.81 
dg =

ne_calc = []
ne_error_calc = []


# Partial derivatives

# ∂(ne)/∂d
partial_d = -(9 * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v0**3) / (g * (p - sigma)))

# ∂(ne)/∂V0
partial_V0 = (9 * d * np.pi * np.sqrt(2) / V0**2) * np.sqrt((n**3 * v0**3) / (g * (p - sigma)))

# ∂(ne)/∂η
partial_n = -(27/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n * v0**3) / (g * (p - sigma)))

# ∂(ne)/∂v0
partial_v0 = -(27/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v0) / (g * (p - sigma)))

# ∂(ne)/∂g
partial_g = (9/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v0**3) / (g**3 * (p - sigma)))

# ∂(ne)/∂ρ
partial_p = (9/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v0**3) / (g * (p - sigma)**3))

# ∂(ne)/∂σ
partial_sigma = -(9/2 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((n**3 * v0**3) / (g * (p - sigma)**3))




for V0, v0 in zip(V0_exp, exp):
    
    ne = - (9 * d * np.pi * np.sqrt(2) / V0) * np.sqrt((eta**3 * v0**3) / (g * (p - sigma)))
    
    delta_ne = np.sqrt(
    (partial_V0 * dV0)**2 +
    (partial_d * dd)**2 +
    (partial_n * dn)**2 +
    (partial_v0 * dv0)**2 +
    (partial_g * dg)**2 +
    (partial_p * dp)**2 +
    (partial_sigma * dsigma)**2)
    
    ne_calc.append(ne)
    ne_error_calc.append(delta_ne)
    

# Results
results = np.column_stack((ne_calc, ne_error_calc))

