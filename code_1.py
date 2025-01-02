import scipy
from scipy import constants
import math

#input
mass_vehicle=2000 #mass of vehicle in kilogram
alpha=2 #slip angle in degree
my=0.2 #friction coefficient

#definitions
camber=2 #camber angle
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

a1_Fx=-21.3
a2_Fx=1144.0
a3_Fx=49.6
a4_Fx=226
a5_Fx=0.069
a6_Fx=-0.006
a7_Fx=0.056
a8_Fx=0.486
a9_Fx=0.028
a10_Fx=0.0
a11_Fx=14.8
a12_Fx=0.022
a13_Fx=0.0

#General Calculations
Fz=scipy.constants.g*mass_vehicle/4000 #load on tire in kilogramm

#Side Force
D_sf=a1_Fy*Fz**2+a2_Fy*Fz #peak factor
C_sf=1.3
B_sf=(a3_Fy*math.sin(a4_Fy*math.atan(a5_Fy*Fz))/C_sf*D_sf)*(1-a12_Fy*abs(camber))
E_sf=a6_Fy*Fz**2+a7_Fy*Fz+a8_Fy
delta_Sh=a9_Fy*camber
delta_Sv=(a10_Fy*Fz**2+a11_Fy*Fz)*camber
phi_sf=(1-E_sf)*(alpha+delta_Sh)+(E_sf/B_sf)*math.atan(B_sf*(alpha+delta_Sh))
Fy=D_sf*math.sin(C_sf*math.atan(B_sf*phi_sf))+delta_Sv

print ("Fy=",Fy)

#Brake Force
D_bf=a1_Fx*Fz**2+a2_Fx*Fz #peak factor
C_bf=1.65
B_bf=(a3_Fx*Fz**2+a4_Fx*Fz)/(C_bf*D_bf*math.e**(a5_Fx*Fz))
E_bf=a6_Fx*Fz**2+a7_Fx*Fz+a8_Fx
phi_bf=(1-E_bf)*alpha+(E_bf/B_bf)*math.atan(B_bf*alpha)
Fx=D_bf*math.sin(C_bf*math.atan(B_bf*phi_bf))
print ("Fx=",Fx)