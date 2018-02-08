import pandas as pd 
import matplotlib.pyplot as plt 
from scipy import interpolate
from scipy import optimize
import numpy as np

import gasdynamics as gd 


class CEA_constants():
	def __init__(self,altitude):
		"""Class that should eventually queery NASA CEA data, for now holds all values constant. Updated 29 Jan 2018.

		Args:
			altitude (float): altitude in meters

		Returns:
			Stores NASA CEA data in atributes
		"""	
		self.gamma = 1.237 #np.mean([1.2534,1.2852])
		self.T_c = 2831.47 # combustion chamber temperature
		self.p_c = 3102640.8 # combustion chamber pressure
		self.rho_c = 3.3826 # combustion chamber density 
		self.a_c = np.sqrt(self.gamma*(1-1/self.gamma)*200.07*self.T_c) # combustion chamber sound speed
		self.Pr = 0.55645 #average throat to exit Prandtl's number
		self.cp = 1.724 #[KJ/KG-K] average throat to exit constant pressure heat capacity
		self.c = 0.003883468 #[millipoise/K^w] viscocity to temperature coefficient
		self.w = 0.678083301 #viscocity to temperature exponent

def clean_pd_series(x_series,y_series):
	"""Takes in two pd series returns cleaned numpy arrays
	"""
	x = x_series.as_matrix();  y= y_series.as_matrix()

	idx_not_nan = ~np.isnan(x)

	x = x[idx_not_nan]; y = y[idx_not_nan]

	x, idx_unique = np.unique(x,return_index=True)
	y = y[idx_unique]

	return (x,y)

def get_area_coord(upper_coord,lower_coord,n):
	x_upper,y_upper = upper_coord
	x_lower,y_lower = lower_coord	

	upper_tck = interpolate.splrep(x_upper,y_upper)
	lower_tck = interpolate.splrep(x_lower,y_lower)

	lower_x_space = np.linspace(x_lower.min(),x_lower.max(),n)
	upper_x_space = np.linspace(x_upper.min(),x_upper.max(),n)

	y_diff = interpolate.splev(upper_x_space,upper_tck) - interpolate.splev(lower_x_space,lower_tck)
	x_diff = upper_x_space - lower_x_space

	theta = np.arctan(y_diff/x_diff)
	#print(theta)
	## TODO: FINISH FINDING AREA COORDINATES (LINE INTEGRAL SA)
	A = np.pi*np.tan(theta)/np.cos(theta)*(interpolate.splev(upper_x_space,upper_tck)**2-interpolate.splev(lower_x_space,lower_tck)**2)
	A[A>10**10] = (interpolate.splev(upper_x_space[A>10**10],upper_tck)**2-interpolate.splev(lower_x_space[A>10**10],lower_tck)**2)*np.pi

	# fig, (ax1,ax2) = plt.subplots(1,2)

	# ax1.plot(x_lower,y_lower,'b-')
	# ax1.plot(x_upper,y_upper,'r-')
	# ax1.axis('equal')
	# ax2.plot(lower_x_space,A)
	# ax2.plot(lower_x_space[A.min()==A],A[A.min()==A],'rx')
	# print(A.min()==A)
	return np.abs(A)
	# plt.axis('equal')

def calc_flow_properties(upper_coord,lower_coord,CEA):
	x_upper,y_upper = upper_coord
	x_lower,y_lower = lower_coord
	# A_t = 
	A = get_area_coord(upper_coord,lower_coord)
	epsilon = A[-1]/A
	def opt_exp(M): gd.expansion_ratio_zero(M,1,CEA.gamma,epsilon) 
	# TODO find zeros of opt_exp to get M
	gd.isentropi

def arc_length_coord(x,y):
	y_dummy = y[1:] - y[:-1]
	x_dummy = x[1:] - x[:-1]
	s_dummy = np.sqrt(y_dummy**2+x_dummy**2)

	s_dummy = np.concatenate((np.array([0]),s_dummy))

	s = np.zeros(s_dummy.shape)

	for i in range(1,s_dummy.shape[0]):
		s[i] = s_dummy[i] + s[i-1]

	return s

if __name__=="__main__":
	n = 1000
	conv_df = pd.read_excel('converge_geom.xlsx')

	# cleaning data...
	x_lower,y_lower = clean_pd_series(conv_df['x (mm)'],conv_df['Copper nozzle y (mm)'])
	x_upper,y_upper = clean_pd_series(conv_df['x (mm)'],conv_df['graphite chamber y (mm)'])

	x_lower = x_lower/1000; y_lower = y_lower/1000
	x_upper = x_upper/1000; y_upper = y_upper/1000

	#interpolating upper curve
	upper_tck = interpolate.splrep(x_upper,y_upper)
	lower_tck = interpolate.splrep(x_lower,y_lower)

	lower_x_space = np.linspace(x_lower.min(),x_lower.max(),n)
	upper_x_space = np.linspace(x_upper.min(),x_upper.max(),n)

	lower_y_space = interpolate.splev(lower_x_space,lower_tck)
	upper_y_space = interpolate.splev(upper_x_space,upper_tck)

	# cea vals
	CEA = CEA_constants(0)
	#expansion_ratio_zero(M_1,M_2,CEA.gamma,epsilon)

	A = get_area_coord((x_upper,y_upper),(x_lower,y_lower),n)

	inv_eps = A/A[-1]
	#print(inv_eps)

	Mach = np.zeros(inv_eps.shape)


	for i in range(len(inv_eps)):

		def M_find(M):
			return inv_eps[i] - 1/M*(2/(CEA.gamma+1)*(1+(CEA.gamma-1)/2*M**2))**((CEA.gamma+1)/(2*(CEA.gamma-1)))
		try:
			Mach[i] = optimize.brentq(M_find,1*10**-7,1)
		except:
			print('except!')
			Mach[i]=0
	print(Mach)
	(T_ratio,p_ratio,rho_ratio,a_ratio) = gd.isentropic_ratios(0,Mach,CEA.gamma)

	#print(inv_eps)
	T = T_ratio*CEA.T_c
	p = p_ratio*CEA.p_c
	rho = rho_ratio*CEA.rho_c
	a = a_ratio*CEA.a_c 

	V = Mach*a

	s = arc_length_coord(lower_x_space,lower_y_space)

	print(type(s))
	csv_array = np.array([lower_x_space,lower_y_space,s,p,T,Mach,A,a,V,rho])
	np.savetxt('converge_contour.csv', csv_array.T, delimiter = ',')
	with open('converge_contour.csv','r') as original: data = original.read()
	with open('converge_contour.csv','w') as modified: modified.write('x,y,s,p,T,M,A,a,V,rho\n' + data)


	plt.plot(lower_x_space,V,'o')

	plt.show()

