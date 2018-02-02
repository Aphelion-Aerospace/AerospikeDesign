import numpy as np 
import matplotlib.pyplot as plt
import scipy.optimize
import gasdynamics as gd
from heat_flux import heat_flux
from plug_nozzle_angelino import plug_nozzle 
import MOC


## NASA CEA CONSTANTS
class CEA_constants():
	def __init__(self,gamma,T_c,p_c,rho_c,a_c,Pr,cp,c,w):
		self.gamma = gamma
		self.T_c = T_c
		self.p_c = p_c
		self.rho_c = rho_c
		self.a_c = a_c
		self.Pr = Pr 
		self.cp = cp 
		self.c = c
		self.w = w

def COST_FNC(design_alt,truncate_ratio,T_w,CEA,r_e,alpha,beta,n):

	### DESIGNING NOZZLE
	(p_atm,T_atm,rho_atm) = gd.standard_atmosphere([design_alt])

	PR = CEA.p_c/p_atm

	M_e = gd.PR_expansion_mach(PR,CEA.gamma)

	expansion_ratio = gd.expansion_ratio(1,M_e,CEA.gamma)#6.64 #8.1273
	# print('Exp. ratio: ' + str(expansion_ratio))
	# print('PR: ' + str(PR))

	A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)

	spike = plug_nozzle(expansion_ratio,A_t,r_e,CEA.gamma,CEA.T_c,CEA.p_c,CEA.a_c,CEA.rho_c,n,truncate_ratio = truncate_ratio)


	### CALCULATING COST
	##	thurst estimation over altitude
	alt_range = np.linspace(0,12000,30)
	(p_atm_r,T_atm_r,rho_atm_r) = gd.standard_atmosphere(alt_range)
	#print(CEA.p_c/p_atm_r)
	thrust_range = np.zeros(alt_range.shape)
	for i in range (alt_range.shape[0]):
		if i==10:
			MOC_mesh = MOC.chr_mesh(spike,gamma,alt_range[i],50,downstream_factor=1.2,plot_chr=0)
		else:
			MOC_mesh = MOC.chr_mesh(spike,gamma,alt_range[i],50,downstream_factor=1.2,plot_chr=0)

		thrust_range[i] = MOC_mesh.compute_thrust('nearest',10)

	work = np.trapz(thrust_range,alt_range)
	plt.plot(alt_range,thrust_range,'o')
	plt.show()
	## heat transfer required

	total_heat_flux = heat_flux(CEA.Pr,CEA.cp,CEA.gamma,CEA.c,CEA.w,CEA.T_c,T_w,spike)

	# print('Work*alpha: ' + str(work*alpha))
	# print('Heat flux*beta: ' + str(total_heat_flux*beta))
	return -alpha*work + total_heat_flux*beta
	


## CONSTANTS OF DESIGN FOR AERODYNAMICS
r_e = 0.067/2 #0.034 # likely too large
n = 1000

## NASA CEA CONSTANTS
gamma = 1.2381 #np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474*10**5
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 


#input variables from NASA CEA in metric units:
Pr=0.55645 #average throat to exit Prandtl's number
cp=1.724 #[KJ/KG-K] average throat to exit constant pressure heat capacity
c=0.003883468 #[millipoise/K^w] viscocity to temperature coefficient
w=0.678083301 #viscocity to temperature exponent

## CONSTANTS OF DESIGN FOR HEAT FLUX
#user input variable in metric units:
T_w=600 #[K] desired temperature of nozzle 


## CONSTANTS OF SIM
alpha = 0.07/8 # 0.07/8 : 1 ratio of alpha : beta gives very similar weights
beta = 0
design_alt = 9000
truncate_ratio = 0.2 # bounds on truncate < 0.1425

CEA = CEA_constants(gamma,T_c,p_c,rho_c,a_c,Pr,cp,c,w)

## CONVERTING TO OPTIMIZABLE FUNCTION
cost_lambda = lambda X: COST_FNC(X[0],X[1],T_w,CEA,r_e,alpha,beta,n)

def min_design_alt(X):
	return X[0] - 3000

def max_design_alt(X):
	return -X[0] + 12000

def min_truncate(X):
	return X[1] - 0.2

def max_truncate(X):
	return -X[1] + 1

cons = [{'type':'ineq', 'fun':min_design_alt},{'type':'ineq', 'fun':max_design_alt},{'type':'ineq', 'fun':min_truncate},{'type':'ineq', 'fun':max_truncate}]
print(cost_lambda([design_alt,truncate_ratio]))

# res = scipy.optimize.minimize(cost_lambda,[design_alt,truncate_ratio],constraints = cons)


# minmizer_kwargs = {'constraints':cons}

# res = scipy.optimize.basinhopping(cost_lambda,[design_alt,truncate_ratio],minimizer_kwargs=minmizer_kwargs)
# print(res)