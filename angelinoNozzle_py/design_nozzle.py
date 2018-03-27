from plug_nozzle_angelino import *
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np
#import aerospike_optimzer



r_e = 0.06 #0.034 # likely too large
expansion_ratio = 11.2 #6.64 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2343# np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 


spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000,truncate_ratio = 0.2)

#tck = interpolate.splrep(spike.x,spike.y)




csv_array = np.array([spike.x,spike.y,spike.s,spike.p,spike.T,spike.M,spike.A,spike.a,spike.V,spike.rho])

np.savetxt('aerospike_contour.csv', csv_array.T, delimiter = ',')




#plt.plot(x_outer,y_outer,x_inner,y_inner)

plt.axis('equal')

plt.show()
