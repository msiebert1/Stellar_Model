import numpy as np
import density_model as dm
from constants import CONST_sigma, CONST_G, CONST_a, CONST_c
import opacity as op
from scipy import interpolate

def ep_correction():
	XHe = np.asarray([.1, .5, .9])
	T7 = np.asarray(np.linspace(0., 4., 17))
	#         0  .25   .5  .75  1    1.25  1.5  1.75  2     2.25  2.5  2.75  3    3.25   3.5   3.75  4
	psi = np.asarray([[1.0, 1.0, 1.0, 1.0, 1.0, 1.05, 1.1, 1.15, 1.25, 1.35, 1.4, 1.43, 1.45, 1.45, 1.45, 1.45, 1.45],
		              [1.0, 1.0, 1.0, 1.0, 1.05, 1.15, 1.4, 1.7, 1.75, 1.7, 1.6, 1.55, 1.5, 1.45, 1.45, 1.45, 1.45],
		              [1.0, 1.0, 1.0, 1.05, 1.15, 1.55, 1.85, 1.9, 1.87, 1.8, 1.7, 1.6, 1.55, 1.5, 1.45, 1.45, 1.45]])

	f = interpolate.RectBivariateSpline(XHe, T7, psi, kx = 2, ky = 2)
	return f

def ep_n(f, XH, XHe, Xcno, T, rho):
	f11 = 1.0
	T9 = T/1.e9
	T7 = T/1.e7
	f11 = 1.0
	g11 = 1 + 3.82*T9 + 1.51*T9*T9 + 0.144*T9*T9*T9 - .0114*T9*T9*T9*T9
	g141 = 1 - 2.*T9 + 3.41*T9*T9 - 2.43*T9*T9*T9
	psi = float(f(XHe, T7))
	ep_pp = 2.57e4*psi*f11*g11*rho*XH*XH*(T9**(-2./3.))*np.exp(-3.381/(T9**(1./3.)))
	ep_cno = 8.24e25*g141*Xcno*XHe*rho*(T9**(-2./3.))*np.exp(-15.231*T9**(-1./3.) - (T9/.8)*(T9/.8))
	return ep_pp, ep_cno

def load1(P_c, T_c, rho_c, m, f, X, Y, Xcno, opacs_interp):
	des_log_T = np.log10(T_c)
	T6 = 1e-6*T_c
	des_log_T6 = np.log10(T6)
	des_log_R = np.log10(rho_c) - 3.*des_log_T6
	log_kappa = float(opacs_interp(des_log_T,des_log_R))
	kappa = 10**log_kappa

	ep_pp, ep_cno = ep_n(f, X, Y, Xcno, T_c, rho_c)
	r = (3.*m/(4.*np.pi*rho_c))**(1./3.)
	l = (ep_pp + ep_cno)*m
	P = P_c + (-3.*CONST_G)*(m**(2./3.))*(4.*np.pi*rho_c/3.)**(4./3.)/(8.*np.pi)

	nabla = (3./(16.*np.pi*CONST_a*CONST_c*CONST_G))*(kappa*l*P_c/(m*T_c**4))

	lnT = np.log(T_c) + (-(np.pi/6.)**(1./3.))*CONST_G*nabla*(rho_c**(4./3.))*(m**(2./3.))/P_c 
	T = np.exp(lnT)
	return r, l, P, T

def load2(R, L, X, Y, Z, M_tot, mu):
	r = R
	l = L
	T = (L/(4.*np.pi*R*R*CONST_sigma))**(1./4.)

	opacs, log_T, log_R = op.read_opacities(X, Y, Z)
	f = interpolate.RectBivariateSpline(log_T,log_R,opacs)
	des_log_T = np.log10(T)
	T6 = 1e-6*T
	des_log_T6 = np.log10(T6)

	P_ks = []
	P_gs = []
	rhos = []

	des_log_rho = -8. + 3.*des_log_T6
	rho_min = 10**des_log_rho
	des_log_rho = 1. + 3.*des_log_T6
	rho_max = 10**des_log_rho

	P_min = dm.Ptot(rho_min, des_log_T, mu)
	P_max = dm.Ptot(rho_max, des_log_T, mu)
	test_Ps = np.linspace(P_min, P_max, num = 100000)
	for el in test_Ps:
		rho = dm.density_rad(np.log10(el), des_log_T, mu)
		des_log_R = np.log10(rho) - 3.*des_log_T6
		log_kappa = float(f(des_log_T,des_log_R))
		kappa = 10**log_kappa

		P_k = (2.*CONST_G*M_tot)/(3.*R*R*kappa)
		P_g = dm.Ptot(rho, des_log_T, mu)
		P_ks.append(P_k)
		P_gs.append(P_g)
		rhos.append(rho)

	P_ks = np.asarray(P_ks)
	P_gs = np.asarray(P_gs)
	diffs = abs(P_ks - P_gs)
	min_loc = diffs.tolist().index(min(diffs))
	P = (P_ks[min_loc] + P_gs[min_loc])/2.
	rho_surf = rhos[min_loc]

	return r, l, P, T

def derivs(y, m, mu, f, X, Y, Xcno , opacs_interp):
	rho = dm.density_rad(np.log10(y[2]), np.log10(y[3]), mu)
	des_log_T = np.log10(y[3])
	T6 = 1e-6*y[3]
	des_log_T6 = np.log10(T6)
	des_log_R = np.log10(rho) - 3.*des_log_T6
	log_kappa = float(opacs_interp(des_log_T,des_log_R))
	kappa = 10**log_kappa

	ep_pp, ep_cno = ep_n(f, X, Y, Xcno, y[3], rho)
	nabla = (3./(16.*np.pi*CONST_a*CONST_c*CONST_G))*(kappa*y[1]*y[2]/(m*y[3]**4))

	drdm = 1./(4.*np.pi*y[0]*y[0]*rho)
	dldm = ep_pp + ep_cno
	dPdm = (-CONST_G*m)/(4.*np.pi*(y[0]**4))
	dTdm = (-CONST_G*m*y[3]*nabla)/(4.*y[2]*np.pi*(y[0]**4))
	return [drdm, dldm, dPdm, dTdm]