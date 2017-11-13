import numpy as np 

# VARIABLE NOMENCLATURE:
#	- M: mach number
#	- nu: prandtl-meyer function

# EXPLICIT RELATIONS ##########################################

# Supersonic flow
def mach_angle(M):
	# computes mach angle (M > 1)
	return np.arcsin(1/M)

def prandtl_meyer(M,gamma):
	# computes prandtl-meyer function (M > 1)
	return np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1))) - np.arctan(np.sqrt(M**2-1))

# Isentropic flow
def isentropic_ratios(M_1,M_2,gamma):
	# computes isentropic flow ratio between two points in flow (T_2/T_1, p_2/p_1, rho_2/rho_1)
	T_ratio = (1+(gamma-1)/2*M_1**2)/(1+(gamma-1)/2*M_2**2)
	p_ratio = T_ratio**(gamma/(gamma-1))
	rho_ratio = T_ratio**(1/(gamma-1))	
	a_ratio = np.sqrt(T_ratio)
	return T_ratio,p_ratio,rho_ratio,a_ratio

def expansion_ratio(M_1,M_2,gamma):
	# returns area of expansion ratio A_2/A_1
	return M_1/M_2*((2+(gamma-1)*M_2**2)/(2+(gamma-1)*M_1**2))**((gamma+1)/(2*(gamma-1)))


# IMPLICIT RELATIONS ##########################################
def prandtl_meyer_zero(M,nu,gamma):
	return np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1))) - np.arctan(np.sqrt(M**2-1)) - nu

def expansion_ratio_zero(M_1,M_2,gamma,epsilon):
	return M_1/M_2*((2+(gamma-1)*M_2**2)/(2+(gamma-1)*M_1**2))**((gamma+1)/(2*(gamma-1))) - epsilon