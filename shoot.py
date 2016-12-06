#shoot.py, author: Matthew Siebert
import numpy as np
from scipy.integrate import odeint
import shoot_routines as sr
import density_model as dm
import matplotlib.pyplot as plt
import opacity as op
from scipy import interpolate

def change_init(guess, d, apply_changes, diff_init, mu, m_start, M_tot, m_right, m_left, X, Y, Z, Xcno, f, opacs_interp):
	"""If apply_changes is True, this function determines the full inward and outward solutions for the initial parameters
	   contained in guess. It then returns these two solutions and the discrepancy vector (diff) which is an array of the 
	   differences in each parameter at the fitting point. If apply_changes is False, this function applies the small change
	   d which contains an independent change to one of the intital parameters. It then uses finite difference to see how the
	   discrepancy vector changes with the current parameter of interest. These derivatives are returned."""
	delta = 0.
	if not apply_changes:
		guess = guess + d
		for el in d:
			if el != 0.:
				delta = el

	rho = dm.density_rad(np.log10(guess[2]), np.log10(guess[3]), mu)

	y_c = sr.load1(guess[2], guess[3], rho, m_start, f, X, Y, Xcno, opacs_interp)
	y_surf = sr.load2(guess[0], guess[1], X, Y, Z, M_tot, mu)
	left = odeint(sr.derivs, y_c, m_left, args = (mu, f, X, Y, Xcno, opacs_interp))
	right = odeint(sr.derivs, y_surf, m_right, args = (mu, f, X, Y, Xcno, opacs_interp))

	left_f = np.asarray([left[:,0][-1], left[:,1][-1], left[:,2][-1], left[:,3][-1]])
	right_f = np.asarray([right[:,0][-1], right[:,1][-1], right[:,2][-1], right[:,3][-1]])

	diff = left_f - right_f

	if not apply_changes:
		derivs = (diff - diff_init)/delta
		return derivs
	else:
		return diff, left, right

def converged(diff, scale):
	"""Check solution for convergence. Returns True if each element of discrepancy vector is within some fraction 
	   rel_tol if the initial guesses."""
	rel_tol = 1.e-4
	if (np.abs(diff[0]) < rel_tol*scale[0] and np.abs(diff[1]) < rel_tol*scale[1] and 
		np.abs(diff[2]) < rel_tol*scale[2] and np.abs(diff[3]) < rel_tol*scale[3]):
		return True
	else:
		return False

def shoot_and_fit(P_c, T_c, R, L, X, Y, Z, m_start, M_tot, m_left, m_right):
	"""Main controller function. Takes initial guesses, applies offsets and calculates how the discrepancy 
	   vector responds in order to populatate the jacobian (jac). A matrix equation is solved in an attempt to 
	   zero the discrepancy vector (diff_init). A fraction of the determines offsets are then used in the next 
	   iteration. This process is repeated until the convergence criteria are met."""
	f = sr.ep_correction()
	opacs, log_T, log_R = op.read_opacities(X, Y, Z)
	opacs_interp = interpolate.RectBivariateSpline(log_T,log_R,opacs)

	mu = dm.mean_molecular_weight(X,Y)
	Xcno = .7*Z

	max_iter = 1000

	guess = np.asarray([R, L, P_c, T_c])
	d = [1.e-2*R, 1.e-2*L, 1.e-2*P_c, 1.e-2*T_c]
	for i in range(max_iter):
		diff_init, left, right = change_init(guess, d, True, None, mu, m_start, M_tot, m_right, m_left, X, Y, Z, Xcno, f, opacs_interp)
		print i, diff_init
		if converged(diff_init, [R, L, P_c, T_c]):
			return left, right, guess 

		dR = [d[0], 0., 0., 0.]
		derivs_1 = change_init(guess, dR, False, diff_init, mu, m_start, M_tot, m_right, m_left, X, Y, Z, Xcno, f, opacs_interp)
		dL = [0., d[1], 0., 0.]
		derivs_2 = change_init(guess, dL, False, diff_init, mu, m_start, M_tot, m_right, m_left, X, Y, Z, Xcno, f, opacs_interp)
		dP = [0., 0., d[2], 0.]
		derivs_3 = change_init(guess, dP, False, diff_init, mu, m_start, M_tot, m_right, m_left, X, Y, Z, Xcno, f, opacs_interp)
		dT = [0., 0., 0., d[3]]
		derivs_4 = change_init(guess, dT, False, diff_init, mu, m_start, M_tot, m_right, m_left, X, Y, Z, Xcno, f, opacs_interp)

		jac = np.transpose([derivs_1, derivs_2, derivs_3, derivs_4])
		d = np.linalg.solve(jac, -1*diff_init)
		d = .2*d
		guess = guess + d

	return "Reached maximum number of iterations"

##Here are the various initial 

# 2 solar masses
# P_c = 1.62e17
# T_c = 2.109e7
# L = 6.998e34
# R = 1.03e11
# X = .7
# Y = .28
# Z = 1 - X - Y
# Xcno = .7*Z
# M_sun = 2.e33
# M_tot = 2.*M_sun

# 3 solar masses
# P_c = 1.148e17
# T_c = 2.347e7
# L = 3.4196e35
# R = 1.276e11
# X = .7
# Y = .28
# Z = 1 - X - Y
# Xcno = .7*Z
# M_sun = 2.e33
# M_tot = 3.*M_sun

# 1.5 solar masses
# P_c = 1.905e17
# T_c = 1.905e7
# L = 2.198e34
# R = 9.151e10
# X = .7
# Y = .28
# Z = 1 - X - Y
# Xcno = .7*Z
# M_sun = 2.e33
# M_tot = 1.5*M_sun

# 10 solar mass
# P_c = 3.715e16
# T_c = 3.048e7
# L = 2.264e37
# R = 2.594e11
# X = .7
# Y = .28
# Z = 1 - X - Y
# Xcno = .7*Z
# M_sun = 2.e33
# M_tot = 10.*M_sun

#30 solar mass
# P_c = 1.9498e16
# T_c = 3.628e7
# L = 4.456e38
# R = 4.853e11
# X = .7
# Y = .28
# Z = 1 - X - Y
# Xcno = .7*Z
# M_sun = 2.e33
# M_tot = 30.*M_sun

#20 Solar Mass
P_c = 2.344e16
T_c = 3.427e7
L = 1.637e38
R = 3.873e11
X = .7
Y = .28
Z = 1 - X - Y
Xcno = .7*Z
M_sun = 2.e33
M_tot = 20.*M_sun

m_start = 1e-9*M_tot
m_left = np.linspace(m_start, .5*M_tot, 1000000)
m_right = np.linspace(M_tot, .5*M_tot, 1000000)
sol_left, sol_right, guess = shoot_and_fit(P_c, T_c, R, L, X, Y, Z, m_start, M_tot, m_left, m_right)
m_left = m_left/M_tot
m_right = m_right/M_tot
print "Done"

#code for calculating and displaying temperature gradient
# lnT = np.log(sol_left[:, 3])
# lnP = np.log(sol_left[:, 2])
# nabla = []
# for i in range(len(m_left)-1):
# 	nabla.append((lnT[i+1] - lnT[i])/(lnP[i+1] - lnP[i]))

# plt.plot(m_left[0:-1], nabla)

# lnT = np.log(sol_right[:, 3])
# lnP = np.log(sol_right[:, 2])
# nabla = []
# for i in range(len(m_right)-1):
# 	nabla.append((lnT[i+1] - lnT[i])/(lnP[i+1] - lnP[i]))

# plt.plot(m_right[0:-1], nabla)
# plt.xlabel('m/M_tot (g)')
# plt.ylabel('dlnT/dlnP')
# plt.savefig('nabla.png')
# plt.show()


print 'R = ', guess[0]
print 'L = ', guess[1]
print 'P_c = ', guess[2]
print 'T_c = ', guess[3]

plt.xlabel('m/M_tot (g)')
plt.ylabel('r (cm)')
plt.plot(m_left, sol_left[:, 0])
plt.plot(m_right, sol_right[:, 0])
plt.xlim([0, 1.05])
# plt.savefig('radius_10.png')
plt.show()

plt.xlabel('m/M_tot (g)')
plt.ylabel('l (erg/s)')
plt.plot(m_left, sol_left[:, 1])
plt.plot(m_right, sol_right[:, 1])
plt.xlim([0, 1.05])
# plt.savefig('luminosity_10.png')
plt.show()

plt.xlabel('m/M_tot (g)')
plt.ylabel('log(P) (g cm^-1 s^-2)')
plt.plot(m_left, np.log10(sol_left[:, 2]))
plt.plot(m_right, np.log10(sol_right[:, 2]))
plt.xlim([0, 1.05])
# plt.savefig('pressure__10.png')
plt.show()

plt.xlabel('m/M_tot (g)')
plt.ylabel('log(T) (K)')
plt.plot(m_left, np.log10(sol_left[:, 3]))
plt.plot(m_right, np.log10(sol_right[:, 3]))
plt.xlim([0, 1.05])
# plt.savefig('temperature_10.png')
plt.show()

