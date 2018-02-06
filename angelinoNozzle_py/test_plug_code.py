from plug_nozzle_angelino import plug_nozzle
import matplotlib.pyplot as plt
import numpy as np
import aerospike_optimizer as ao 

r_e = 0.027 
expansion_ratio = 6.64 #8.1273, based on pressure ratio
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2343# np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c)
T_w = 600
alpha = 1
beta = 1
truncate_ratio_init = 0.2
design_alt_init = 9144 # 30 % greater

opt_aero = ao.aerospike_optimizer(r_e,T_w,alpha,beta,design_alt_init,truncate_ratio_init,chr_mesh_n=120,no_alt_range = 30,no_core=1)

#design diverging section, 20% truncation

plug1 = opt_aero.spike_init

plug1.define_compression(1.15/1000,4.51/1000,1,12.91/1000,10000)

plug1.plot_contour(plt)
plt.axis('equal')
plt.show()

plug1.save_to_csv()
