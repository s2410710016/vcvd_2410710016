#%%

import scipy
from scipy import constants
import numpy as np
import matplotlib.pyplot as plt

#input
mass_vehicle=2000 #mass of vehicle in kilogram

#definitions
camber=0 #camber angle
a1_Fy=-22.1
a2_Fy=1011.0
a3_Fy=1078.0
a4_Fy=1.82
a5_Fy=0.208
a6_Fy=0.0
a7_Fy=-0.354
a8_Fy=0.707
a9_Fy=0.028
a10_Fy=0.0
a11_Fy=14.8
a12_Fy=0.022
a13_Fy=0.0


#General Calculations
Fz=8
#scipy.constants.g*mass_vehicle/4000 #load on each tire in kilogramm

print ("Fz=",Fz)

#Side Force
D_sf=a1_Fy*Fz**2+a2_Fy*Fz #peak factor
C_sf=1.3
B_sf=(a3_Fy*np.sin(a4_Fy*np.arctan(a5_Fy*Fz))/(C_sf*D_sf))*(1-a12_Fy*np.abs(camber))
E_sf=a6_Fy*Fz**2+a7_Fy*Fz+a8_Fy
delta_Sh=a9_Fy*camber
delta_Sv=(a10_Fy*Fz**2+a11_Fy*Fz)*camber

def phi_sf(alpha):
    return(1-E_sf)*(alpha+delta_Sh)+(E_sf/B_sf)*np.arctan(B_sf*(alpha+delta_Sh))

def Fy0(alpha):
    return D_sf*np.sin(C_sf*np.arctan(B_sf*phi_sf(alpha)))+delta_Sv

alpha=np.linspace(-25,25,100000)

print ("D",D_sf)
print ("E",E_sf)
print ("C", C_sf)
print ("B", B_sf)

print (Fy0(1))
print (Fy0(15))

plt.plot(alpha, phi_sf(alpha))
plt.plot(alpha, Fy0(alpha))
plt.show()
