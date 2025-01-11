import scipy
from scipy import constants
import numpy as np
import matplotlib.pyplot as plt

def f(alpha):
    return 3*np.arctan(alpha)

alpha=np.linspace(-10,10,1000)

plt.plot(alpha,f(alpha))
plt.show()