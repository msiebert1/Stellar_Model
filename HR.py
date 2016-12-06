#HR.py, author: Matthew Siebert
import numpy as np
import matplotlib.pyplot as plt
from constants import CONST_sigma

#Creates an HR diagram for the stars modeled in the paper
L = np.asarray([1.84e34, 6.10e34, 2.98e35, 2.06e37, 1.46e38, 3.79e38])
R = np.asarray([1.03e11, 1.13e11, 1.38e11, 2.69e11, 3.82e11, 4.55e11])
T = (L/(4.*np.pi*R*R*CONST_sigma))**(1./4.)
M = np.asarray([1.5, 2., 3., 10., 20., 30])
M = 2.e33*M

logT = np.log10(T)
logL = np.log10(L/3.9e33)
logM = np.log10(M/2.e33)

p = np.polyfit(logT,logL, 1)
x = np.linspace(3.5, 5., 100000)
y = p[0]*x + p[1]
print p
print T 
print np.log10(L/3.9e33)
plt.plot(logT, logL, 'bo')
plt.plot(x,y, 'r--')
ax = plt.gca()
ax.invert_xaxis()
plt.xlabel("log(Teff) (K)")
plt.ylabel("log(L) (in L_sun)")
plt.show()
