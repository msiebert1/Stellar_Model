import numpy as np
import density_model as dm
import shoot_routines as sr

T6_1 = 10.
T6_2 = 15.
T6_3 = 25.

T7_1 = T6_1/10.
T7_2 = T6_2/10.
T7_3 = T6_3/10.

def W(Zj, Zk, Aj, Ak):
	return Zj*Zj*Zk*Zk*(Aj*Ak/(Aj +Ak))

def Tau(W, T7):
	return 19.721*W**(1./3.)*T7**(-1./3.)

def E0(W,T7):
	return 5.665*W**(1./3.)*T7**(2./3.)

def nu(W,T7):
	return 6.574*W**(1./3.)*T7**(-1./3.) - 2./3.

print "Problem #1:"
W1 = W(1.,1.,1.,1.)
nu_pp1 = nu(W1, T7_1)
nu_pp2 = nu(W1, T7_2)
E0_pp1 = E0(W1, T7_1)
E0_pp2 = E0(W1, T7_2)
print "Part a:"
print "nu (T6 = 10, 15) = ", nu_pp1, ", ", nu_pp2
print "E0 (T6 = 10, 15) (keV) = ", E0_pp1,", ",  E0_pp2
print


W2 = W(4., 1., 7., 1.)
nu_BeB1 = nu(W2, T7_2)
E0_BeB1 = E0(W2, T7_2)
print "Part b:"
print "nu (T6 = 15) = ", nu_BeB1
print "E0 (T6 = 15) (keV) = ", E0_BeB1
print 

W3 = W(7.,1., 14.,1.)
nu_NO1 = nu(W3, T7_2)
nu_NO2 = nu(W3, T7_3)
E0_NO1 = E0(W3, T7_2)
E0_NO2 = E0(W3, T7_3)
print "Part c:"
print "nu (T6 = 15, 25) = ", nu_NO1, ", ", nu_NO2
print "E0 (T6 = 15, 25) (keV) = ", E0_NO1,", ", E0_NO2
print
print

print "Problem #2:"
mu_exterior = dm.mean_molecular_weight(1.,0.)
mu_interior = dm.mean_molecular_weight(0.,1.)
print "For an ideal gas mu is proprotional to density. \nTherefore, rho_int/rho_ext = mu_int/rho_ext"
print "rho_int/rho_ext = ", mu_interior/mu_exterior
print 

print "Problem #3:"
X = .7
Y = .28
Z = 1 - X - Y
Xcno = .7*Z
T = 18*1.e6
rho = 80

f = sr.ep_correction()
ep_pps = []
ep_cnos = []
ep_pp, ep_cno = sr.ep_n(f, X, Y, Xcno, T, rho)
print "epsillon_pp/epsilon_cno = ", ep_pp/ep_cno
