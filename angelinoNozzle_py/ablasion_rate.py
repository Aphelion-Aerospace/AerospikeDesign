import numpy as np 
import gasdynamics as gd 
from scipy import optimize
from angelino_nozzle_design import plug_nozzle
import matplotlib.pyplot as plt

def ideal_thrust(A_t,A_e,gamma,p_c,T_c,rho_c,a_c,altitude,m_dot):

	p_atm,T_atm,rho_atm = gd.standard_atmosphere([altitude])

	epsilon = A_e/A_t
	
	M_e = optimize.fsolve(lambda M: gd.expansion_ratio_zero(1,M,gamma,epsilon),5)
	# print(M_e)
	# throat conditions
	T_ratio,p_ratio,rho_ratio,a_ratio = gd.isentropic_ratios(0,1,gamma)
	T_t = T_ratio*T_c; p_t = p_ratio*p_c; rho_t = rho_ratio*rho_c; a_t = a_ratio*a_c;
	
	# exit conditions
	T_ratio,p_ratio,rho_ratio,a_ratio = gd.isentropic_ratios(0,M_e,gamma)
	T_e = T_ratio*T_c; p_e = p_ratio*p_c; rho_e = rho_ratio*rho_c; a_e = a_ratio*a_c;

	# mass flow rate (constant in ideal case)
	# m_dot = rho_t*A_t*a_t*1 
	# print(m_dot)

	# print('rho calc ' + str(rho_t))
	# print('a_t calc ' + str(a_t))
	# print('V_e calc ' + str(M_e*a_e))
	# print('p_e - p_atm' + str(p_e - p_atm))
	# print('M_e ' + str(M_e))

	return m_dot*M_e*a_e #+ (p_e - p_atm)*A_e

def A_t(m,x_a,y_a,x_b):
	S = 2*np.pi*np.sqrt(1 + m**2) * (m/2*(x_b**2-x_a**2) + (y_a - m*x_a)*(x_b-x_a))	

	return S 

def S_o_t(m,x_a,y_a,x_b,alpha,t):
	theta = np.arctan(m)

	xa_t = x_a -alpha*t*np.cos(theta)
	ya_t = y_a -alpha*t*np.sin(theta)
	xb_t = x_b +alpha*t*np.cos(theta)

	return A_t(m,xa_t,ya_t,xb_t)

def ideal_thrust_vector(alpha,t_vec,plug1,A_e,gamma,p_c,T_c,rho_c,a_c,altitude):

	m_dot = plug1.rho[0]*plug1.A_t*plug1.V[0]
	thrust_vec = []
	for t in t_vec:
		A_t_alpha = S_o_t(m,plug1.x[0],plug1.y[0],plug1.lip_x,alpha,t)
		thrust_alpha = ideal_thrust(A_t_alpha,A_e,gamma,p_c,T_c,rho_c,a_c,altitude,m_dot)
		thrust_vec.append(thrust_alpha)

	return thrust_vec

def num_runs(thrust_vector,cutoff):
	initial = thrust_vector[0]
	final = thrust_vector[-1]
	
	percent_drop = (initial - final)/initial

	print(percent_drop)
	num = np.floor(cutoff/percent_drop)

	return int(num)

if __name__ == '__main__':


	# engine constants
	r_e = 0.027
	A_e = np.pi*r_e**2
	gamma = 1.237 #np.mean([1.2534,1.2852])
	T_c = 2831.47 # combustion chamber temperature
	p_c = 3102640.8 # combustion chamber pressure
	rho_c = 3.3826 # combustion chamber density 
	a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) # combustion chamber sound speed

	altitude = 9144

	p_atm,T_atm,rho_atm = gd.standard_atmosphere([altitude])

	plug1 = plug_nozzle(altitude,r_e,gamma,T_c,p_c,a_c,rho_c,100,truncate_ratio = 1)

	m = (plug1.lip_y - plug1.y[0])/(plug1.lip_x - plug1.x[0])

	# simulation constant 
	num_t = 10

	t_vec = np.linspace(0,4,num_t)
	
	alpha = 0.03/1000
	beta = 0.05/1000


	thrust_vector_alpha = ideal_thrust_vector(alpha,t_vec,plug1,A_e,gamma,p_c,T_c,rho_c,a_c,altitude)
	thrust_vector_beta = ideal_thrust_vector(beta,t_vec,plug1,A_e,gamma,p_c,T_c,rho_c,a_c,altitude)
	
	thrust_vector_beta[-1]
	fig, ax1 = plt.subplots(1,1)

	alpha_runs = num_runs(thrust_vector_alpha,0.1)
	beta_runs = num_runs(thrust_vector_beta,0.1)

	print('Max runs: ' + str(alpha_runs))
	print('Min runs: ' + str(beta_runs))

	ax1.plot(t_vec,thrust_vector_alpha,label = "Best case")
	ax1.plot(t_vec,thrust_vector_beta, label = "Worst case")
	ax1.set_xlabel('Burn time (s)')
	ax1.set_ylabel('Thrust performance (N)')
	ax1.set_title('Ablasion Peformance for Single Fire')
	plt.legend()
	fig.set_size_inches(18.5,10.5)
	plt.savefig('final_report/ablasion',dpi=100) 
	# print(thrust_vector)
	plt.show()


	# A_t = S_o_t(m,plug1.x[0],plug1.y[0],plug1.lip_x,alpha,0)

	# throat_areas = np.linspace(A_t/10,A_t*10,10)

	# thrust = ideal_thrust(throat_areas,A_e,gamma,p_c,T_c,rho_c,a_c,altitude)
	# # thrust2 = plug1.calc_ideal_thrust(p_atm)
	# print(thrust)
	# print(thrust2)
	## TODO: Calculate throat area based on time since fire