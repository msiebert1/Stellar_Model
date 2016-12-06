#density_model.py, author: Matthew Siebert
import numpy as np
from constants import CONST_a, CONST_Rgas


def mean_molecular_weight(X,Y):
    """Returns the mean molecular weight of a fully ionized gas given and 
       composition (X and Y)."""
    return 2./(1 + 3.*X + Y/2.)

def density(logP,logT,mu):
	"""Returns the density of an ideal gas given logP, logT, and 
	   mean molecular weight mu"""
	P = 10.**logP
	T = 10.**logT
	return (mu/CONST_Rgas)*(P/T)

def density_rad(logP,logT,mu):
	"""Returns the density of an ideal gas accounting for radiation pressure 
	   given logP, logT, and mean molecular weight mu."""
	P = 10.**logP
	T = 10.**logT
	return (mu/CONST_Rgas)*((P/T) - (CONST_a*T**3)/3.)

def Pgas(rho,logT,mu):
	"""Returns the pressure of an ideal gas given density rho, logT, and 
	   mean molecular weight mu"""
	return (CONST_Rgas/mu)*rho*(10**(logT))

def Ptot(rho, logT, mu):
	"""Returns the total pressure of an ideal gas accounting radiation 
	   pressure given density rho, logT, and mean molecular weight mu"""
	T = 10.**logT
	return Pgas(rho,logT,mu) + (CONST_a/3.)*(T**4.)

# print 'Part (a):'
# mu = mean_molecular_weight(0.,.98)
# rho = density_rad(16.85, 7.55, mu)
# print 'density = ', rho, ' g*cm^-3'
# beta = Pgas(rho,7.55,mu)/(10**(16.85))
# print 'beta = ', beta

# print

# print 'Part (b): '
# mu = mean_molecular_weight(.7, .28)
# rho = density_rad(16.87, 6.91, mu)
# print 'density = ', rho, ' g*cm^-3'
# beta = Pgas(rho,6.91,mu)/(10**(16.87))
# print 'beta = ', beta