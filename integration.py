import numpy as np
from scipy.integrate import odeint
import shoot_routines as sr
import density_model as dm
import matplotlib.pyplot as plt
import opacity as op
from scipy import interpolate

P_c = 2.7e16
T_c = 1.6e7
L = 3.9e33
R = 7.e10

X = .7
Y = .28
Z = 1 - X - Y
Xcno = .7*Z
M_tot = 2.e33
nabla = .4

f = sr.ep_correction()
opacs, log_T, log_R = op.read_opacities(X, Y, Z)
opacs_interp = interpolate.RectBivariateSpline(log_T,log_R,opacs)

mu = dm.mean_molecular_weight(X,Y)
rho_c = dm.density_rad(np.log10(P_c), np.log10(T_c), mu)

m_start = 1e-9*M_tot
y_c = sr.load1(P_c, T_c, rho_c, m_start, f, X, Y, Xcno, opacs_interp)
y_surf = sr.load2(R, L, X, Y, Z, M_tot, mu)

m_left = np.linspace(m_start, M_tot/2., 1000000)
sol_left = odeint(sr.derivs, y_c, m_left, args = (mu, f, X, Y, Xcno, opacs_interp))
m_right = np.linspace(M_tot, M_tot/2., 1000000)
sol_right = odeint(sr.derivs, y_surf, m_right, args = (mu, f, X, Y, Xcno, opacs_interp))

print "Fitting point values (cgs):"
print "From center:"
print "r = ", sol_left[:,0][-1], "l = ", sol_left[:,1][-1], "P = ", sol_left[:,2][-1], "T = ", sol_left[:,3][-1]
print "From surface:"
print "r = ", sol_right[:,0][-1], "l = ", sol_right[:,1][-1], "P = ", sol_right[:,2][-1], "T = ", sol_right[:,3][-1]

plt.xlabel('m (g)')
plt.ylabel('r (cm)')
plt.plot(m_left, sol_left[:, 0])
plt.plot(m_right, sol_right[:, 0])
plt.show()

plt.xlabel('m (g)')
plt.ylabel('l (erg/s)')
plt.plot(m_left, sol_left[:, 1])
plt.plot(m_right, sol_right[:, 1])
plt.show()

plt.xlabel('m (g)')
plt.ylabel('log(P) (g cm^-1 s^-2)')
plt.plot(m_left, np.log10(sol_left[:, 2]))
plt.plot(m_right, np.log10(sol_right[:, 2]))
plt.show()

plt.xlabel('m (g)')
plt.ylabel('log(T) (K)')
plt.plot(m_left, np.log10(sol_left[:, 3]))
plt.plot(m_right, np.log10(sol_right[:, 3]))
plt.show()