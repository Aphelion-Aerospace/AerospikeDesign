import numpy as np 
import gasdynamics as gd
#from heat_flux import heat_flux
from plug_nozzle_angelino import plug_nozzle 
import MOC



## NASA CEA CONSTANTS
class CEA_constants():
	def __init__(self,gamma,T_c,p_c,rho_c,a_c):
		self.gamma = gamma
		self.T_c = T_c
		self.p_c = p_c
		self.rho_c = rho_c
		self.a_c = a_c


def COST_FNC(design_alt,truncate_ratio,CEA,r_e,alpha,beta,n):
	(p_atm,T_atm,rho_atm) = gd.standard_atmosphere(design_alt)

	PR = CEA.p_c/p_atm

	M_e = gd.PR_expansion_mach(PR,CEA.gamma)

	expansion_ratio = gd.expansion_ratio(1,M_e,CEA.gamma)#6.64 #8.1273

	A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)

	#### DESIGN OF SPIKE

	spike = plug_nozzle(expansion_ratio,A_t,r_e,CEA.gamma,CEA.T_c,CEA.p_c,CEA.a_c,CEA.rho_c,n,truncate_ratio = truncate_ratio)

	# flip spike!!
	spike.y = spike.y*-1
	spike.lip_y = spike.lip_y*-1


	MOC_mesh = MOC.chr_mesh(spike,gamma,12000,50,downstream_factor=1.2,plot_chr=1)

	thrust = MOC_mesh.compute_thrust('nearest',10)

	return alpha*thrust 
# def COST(...alpha,beta):
# 	# spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 1)
# 	# heat = heat_flux(Pr,Cp,Gamma,c,w,To,spike)
# 	# thrust = MOC.thrust(...)
# 	return alpha*heat + beta*thurst
# 	# Pr = ??, c = ??, w = ??, To = Tc??

## CONSTANTS OF DESIGN
r_e = 0.067/2 #0.034 # likely too large
n = 1000

## NASA CEA CONSTANTS
gamma = 1.2381 #np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474*10**5
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 

## CONSTANTS OF SIM
alpha = 1
beta = 0.2
design_alt = 9000
truncate_ratio = 0.2

CEA = CEA_constants(gamma,T_c,p_c,rho_c,a_c)

print(COST_FNC(design_alt,truncate_ratio,CEA,r_e,alpha,beta,n))