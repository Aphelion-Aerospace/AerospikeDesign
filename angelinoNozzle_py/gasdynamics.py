import numpy as np 
from copy import copy
#### GASDYNAMICS V1.0

# VARIABLE NOMENCLATURE:
#	- M: mach number
#	- nu: prandtl-meyer function
#	- T: temperature
#	- p: pressure
#	- rho: density

# EXPLICIT RELATIONS ##########################################

# Supersonic flow
def mach_angle(M,degrees=0):
	# computes mach angle (M > 1)
	mu = np.arcsin(1/M)
	if (degrees):
		return mu*180/np.pi
	else:
		return mu

def prandtl_meyer(M,gamma,degrees=0):
	# computes prandtl-meyer function (M > 1)
	nu = np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1))) - np.arctan(np.sqrt(M**2-1))
	if (degrees):
		return nu*180/np.pi 
	else:
		return nu

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

# Important rocket relations
def PR_expansion_mach(PR,gamma):
	# returns mach number given pressure ratio
	return np.sqrt(((PR)**((gamma-1)/gamma)-1)*2/(gamma-1))

# IMPLICIT RELATIONS ##########################################
def prandtl_meyer_zero(M,nu,gamma):
	return np.sqrt((gamma+1)/(gamma-1))*np.arctan(np.sqrt((gamma-1)/(gamma+1)*(M**2-1))) - np.arctan(np.sqrt(M**2-1)) - nu

def expansion_ratio_zero(M_1,M_2,gamma,epsilon):
	return M_1/M_2*((2+(gamma-1)*M_2**2)/(2+(gamma-1)*M_1**2))**((gamma+1)/(2*(gamma-1))) - epsilon

def mach_angle_velocity_ratio(mu,W,gamma):
	return np.sin(mu) - np.sqrt((gamma-1)/2*(1/W**2-1))

# STANDARD ATMOSPHERE #########################################
# constants declaration


def standard_atmosphere(altitude):
	# 1976 US Standard Atmosphere

	ps = 101325
	rhos = 1.225
	Ts = 288.15

	R = 287.0531
	g_0 = 9.80665
	a = -6.5/1000
	earth_radius = 6.356766*10**6
	altitude = copy(altitude)
	p_atm_range = copy(altitude)
	T_atm_range = copy(altitude)
	rho_atm_range = copy(altitude)
	for i in range(len(altitude)):
		h = earth_radius/(earth_radius + altitude[i])*altitude[i] #geopotential alt.
		#h = 20000

		if (h<=11000):
		# gradient (troposphere) region
			T_atm = Ts + a*h;
			p_atm = ps*(T_atm/Ts)**-(g_0/(a*R))
			rho_atm = rhos*(T_atm/Ts)**-(1+g_0/(a*R))
		else:
		# contant temperature (stratosphere)region ~20000m (65616.8ft), should check
			# values at boundary
			T_atm = Ts + a*11000
			p1 = ps*(T_atm/Ts)**(-g_0/(a*R))
			rho1 = rhos*(T_atm/Ts)**-(1+g_0/(a*R))

			# values at elevation
			p_atm = p1*np.exp(-g_0/(R*T_atm)*(h-11000))
			rho_atm = rho1*np.exp(-g_0/(R*T_atm)*(h-11000))
		p_atm_range[i] = p_atm
		T_atm_range[i] = T_atm 
		rho_atm_range[i] = rho_atm 

	if len(altitude) > 1:
		return (p_atm_range,T_atm_range,rho_atm_range)
	else:
		return (p_atm,T_atm,rho_atm)