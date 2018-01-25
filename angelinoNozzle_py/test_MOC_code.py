import numpy as np
import matplotlib.pyplot as plt
from plug_nozzle_angelino import plug_nozzle
from MOC import chr_mesh

r_e = 0.067/2 #0.034 # likely too large
expansion_ratio = 8.67 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2381 #np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474*10**5
rho_c = 3.3826
R = (1-1/gamma)*1.8292*1000#1.8292 1.6196
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 
n = 1000

# design plug nozzle
spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 0.5)

MOC_mesh = chr_mesh(spike,gamma,5278,30,downstream_factor=1.2,plot_chr=1)

spike.plot_contour(plt)
plt.plot(MOC_mesh.x[MOC_mesh.ID_contour_chr],MOC_mesh.y[MOC_mesh.ID_contour_chr],'bx')

# plt.plot(spike.x,spike.rho,'b--',MOC_mesh.x[MOC_mesh.ID_contour_chr],MOC_mesh.rho[MOC_mesh.ID_contour_chr],'r-')

plt.axis('equal')
plt.show()