from plug_nozzle_angelino import plug_nozzle
from MOC import chr_mesh
import numpy as np
import matplotlib.pyplot as plt

n = 50
num_its =30
r_e = 0.067/2 #0.034 # likely too large
expansion_ratio = 6.64 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2343# np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 


spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000,truncate_ratio = 1)

spike.plot_contour(plt)


MOC_mesh = chr_mesh(spike,gamma,9000,n,num_its,downstream_factor=1.1,plot_chr=1)

plt.show()