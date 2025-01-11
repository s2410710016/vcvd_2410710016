

import scipy
from scipy import constants
import math
import numpy as np
import matplotlib.pyplot as plt

alpha_grad=2
alpha=alpha_grad*np.pi/180

#Sigma definieren
def Sigma_x(kappa):
    return -kappa/(1+kappa)

def Sigma_y(kappa):
    return -np.tan(alpha)/(1+kappa)

def Sigma(kappa):
    return np.sqrt(Sigma_x(kappa)**2+Sigma_y(kappa)**2)

kappa=np.linspace(0,100,100)

print(Sigma(0))
print(Sigma_x(0))
print(Sigma_y(0))

print(Sigma(100))
print(Sigma_x(100))
print(Sigma_y(100))

#Diagramm plotten
plt.plot(kappa, Sigma_x(kappa), label='Sigma_x')
plt.plot(kappa, Sigma_y(kappa), label='Sigma_y')
plt.plot(kappa, Sigma(kappa), label='Sigma')
plt.legend(loc='best')
plt.show()

x=np.tan(alpha)
print("tangens alpha", x)
