from plug_nozzle_angelino import plug_nozzle
import matplotlib.pyplot as plt
import numpy as np
r_e = 0.067/2 
expansion_ratio = 6.64 #8.1273, based on pressure ratio
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2343# np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 

#design diverging section, 20% truncation
plug1 = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000,truncate_ratio=0.2)

plug1.define_compression(1.15/1000,4.51/1000,1,12.91/1000,10000)

plug1.plot_contour(plt)
plt.axis('equal')
plt.show()

plug1.save_to_csv()
