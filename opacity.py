import numpy as np
from scipy import interpolate


def read_opacities(myX, myY, myZ):
	f = open('opacity_tables.txt')
	lines = f.readlines()
	tables = lines[62:188]

	#desired composition
	X = myX
	Y = myY
	Z = myZ

	for line in tables:
		if float(line.split()[5][2:]) == X and float(line.split()[6][2:]) == Y:
			table_num = line.split()[2]

	print "Using table ", table_num
	data = lines[241:]

	for i in range(len(data)):
		if len(data[i].split()) > 3 and data[i].split()[2] == table_num:
			loc = i

	opac_data = data[loc+4:loc+76]

	log_R = np.asarray(opac_data[0].split()[1:]).astype(float)
	log_T = []
	for line in opac_data[2:]:
		if len(line.split()[1:]) == 19:
			log_T.append(line.split()[0])

	log_T = np.asarray(log_T).astype(float)

	opacs = []
	for line in opac_data[2:]:
		if len(line.split()[1:]) == 19:
			opacs.append(np.asarray(line.split()[1:]).astype(float))

	opacs = np.asarray(opacs)
	return opacs, log_T, log_R


# opacs, log_T, log_R = read_opacities(.5, .5, 0.)

# des_log_T = 6.3
# T = 10**(des_log_T)
# T6 = 1e-6*T
# des_log_T6 = np.log10(T6)
# des_log_rho = .3
# des_log_R = des_log_rho - (3*des_log_T6)

# f = interpolate.RectBivariateSpline(log_T,log_R,opacs)
# log_kappa = float(f(des_log_T,des_log_R))
# print 'Part (a):'
# print 'kappa = ', 10**log_kappa, ' cm^2*g^-1'


# des_log_T = 5.
# T = 10**(des_log_T)
# T6 = 1e-6*T
# des_log_T6 = np.log10(T6)
# des_log_rho = -4.
# des_log_R = des_log_rho - (3*des_log_T6)

# log_kappa = float(f(des_log_T,des_log_R))
# print
# print 'Part (b):'
# print 'kappa = ', 10**log_kappa, ' cm^2*g^-1'