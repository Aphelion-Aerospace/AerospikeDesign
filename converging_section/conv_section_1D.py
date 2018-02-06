import pandas as pd 
import matplotlib.pyplot as plt 
from scipy import interpolate
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
		self.gamma = 1.2381 #np.mean([1.2534,1.2852])
		self.T_c = 2833.63 # combustion chamber temperature
		self.p_c = 34.474*10**5 # combustion chamber pressure
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

if __name__=="__main__":
	conv_df = pd.read_excel('converge_geom.xlsx')

	# cleaning data...
	x_lower,y_lower = clean_pd_series(conv_df['x'],conv_df['y'])
	x_upper,y_upper = clean_pd_series(conv_df['x.1'],conv_df['y.1'])

	#interpolating upper curve
	upper_tck = interpolate.splrep(x_upper,y_upper)


	# cea vals
	CEA = CEA_constants(altitude)
	expansion_ratio_zero(M_1,M_2,CEA.gamma,epsilon)


	plt.plot(x_lower,y_lower)
	plt.plot(x_upper,y_upper)
	plt.show()

