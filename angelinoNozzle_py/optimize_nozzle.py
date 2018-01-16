import numpy as np 
from heat_flux import heat_flux
from plug_nozzle_angelino import plug_nozzle 
import MOC

def COST(...alpha,beta):
	spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 1)
	heat = heat_flux(Pr,Cp,Gamma,c,w,To,spike)
	thrust = MOC.thrust(...)
	return alpha*heat + beta*thurst
	# Pr = ??, c = ??, w = ??, To = Tc??

