#%%

import scipy
from scipy import constants
import numpy as np
import matplotlib.pyplot as plt
import argparse

#argparse

import argparse

def parse_config(file_path):
    variables = {}
    with open(file_path, "r") as file:
        for line in file:
            key, value = line.strip().split("=")
            variables[key] = eval(value)  # Vorsicht: eval kann unsicher sein!
    return variables

def main():
    parser = argparse.ArgumentParser(description="Lade Variablen aus einer Datei.")
    parser.add_argument("config_file", type=str, help="Pfad zur Konfigurationsdatei (.txt)")
    
    args = parser.parse_args()
    
    config = parse_config(args.config_file)
    
    print("Geladene Variablen:")
    for key, value in config.items():
        print(f"{key}: {value}")

if __name__ == "__main__":
    main()

#input
mass_vehicle=2000 #mass of vehicle in kilogram
alpha_grad=2 #slip angle
alpha=alpha_grad*np.pi/180 #slip angle
my=0.2 #friction coefficient

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

kappa=np.linspace(0,1,100)

#General Calculations
#Fz=scipy.constants.g*mass_vehicle/4000 #load on each tire in kilogramm
Fz=4
print ("Fz=",round(Fz,2),"kN")

#Side Force
D_sf=a1_Fy*Fz**2+a2_Fy*Fz #peak factor
C_sf=1.3 #shape factor
B_sf=(a3_Fy*np.sin(a4_Fy*np.arctan(a5_Fy*Fz))/(C_sf*D_sf))*(1-a12_Fy*np.abs(camber)) #stiffness factor
E_sf=a6_Fy*Fz**2+a7_Fy*Fz+a8_Fy #curvature factor
delta_Sh=a9_Fy*camber #horizontal shift
delta_Sv=(a10_Fy*Fz**2+a11_Fy*Fz)*camber #vertical shift

def phi_sf(alpha):
    return(1-E_sf)*(alpha_grad+delta_Sh)+(E_sf/B_sf)*np.arctan(B_sf*(alpha+delta_Sh))

def Fy(alpha):
    return D_sf*np.sin(C_sf*np.arctan(B_sf*phi_sf(alpha)))+delta_Sv

print ("Fy=",round(Fy(alpha_grad),2),"N bei Alpha: ", alpha_grad,"°")

#Brake Force
D_bf=a1_Fx*Fz**2+a2_Fx*Fz #peak factor
C_bf=1.65 #shape factor
B_bf=(a3_Fx*Fz**2+a4_Fx*Fz)/(C_bf*D_bf*np.e**(a5_Fx*Fz)) #stiffness factor
E_bf=a6_Fx*Fz**2+a7_Fx*Fz+a8_Fx #curvature factor

def phi_bf(kappa):
    return (1-E_bf)*kappa+(E_bf/B_bf)*np.atan(B_bf*kappa)*180/np.pi
def Fx(kappa):
    return D_bf*np.sin(C_bf*np.atan(B_bf*phi_bf(kappa)))
#print ("Fx=",round(Fx(alpha_grad),2),"N bei Alpha: ", alpha_grad, "°")

#Sigma definieren
def Sigma_x(kappa):
    return -kappa/(1+kappa)

print ("Sigma_x:", Sigma_x(1))

def Sigma_y(kappa):
    return -np.tan(alpha)/(1+kappa)

print ("Sigma_y:", Sigma_y(1))

def Sigma(kappa):
    return np.sqrt(Sigma_x(kappa)**2+Sigma_y(kappa)**2)


#Funktionen Fx und Fy definieren
def Fx_new(kappa):
    return -(Sigma_x(kappa)/Sigma(kappa))*Fx(kappa)

def Fy_new(kappa):
    return -(Sigma_y(kappa)/Sigma(kappa))*Fy(alpha)

#Kappa als Laufvariable

#alpha=np.linspace(-20,20,100)

#Diagramm plotten
plt.plot(kappa*100, Fx_new(kappa), label='Fx')
plt.plot(kappa*100, Fy_new(kappa), label='Fy')

#plt.plot(alpha, Fy(alpha), label='Fy_grad')
#plt.plot(alpha, Fx(alpha), label='Fx_grad')
plt.legend(loc='best')
plt.show()

