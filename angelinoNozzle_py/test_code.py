from plug_nozzle_angelino import plug_nozzle
from MOC import chr_mesh
import gasdynamics as gd 
import numpy as np 
import matplotlib.pyplot as plt

#### NOZZLE INITIAL CONDITIONS
r_e = 0.067/2 #0.034 # likely too large
expansion_ratio = 6.64 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2381 #np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474*10**5
rho_c = 3.3826
R = (1-1/gamma)*1.8292*1000#1.8292 1.6196
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 
n = 1000


spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 1.0)
alts = np.linspace(0,15000,100)
p_atms,T_atms,rho_atms = gd.standard_atmosphere(alts)

thrusts = np.zeros(p_atms.shape)

for i in range(len(p_atms)):
	thrusts[i] = spike.calc_ideal_thrust(p_atms[i])

plt.plot(alts,thrusts)


# MOC_mesh = chr_mesh(spike,gamma,0,50,downstream_factor=1.2,plot_chr=1)

#plt.plot(MOC_mesh.x,MOC_mesh.y,'ro')

#plt.plot(MOC_mesh.x[MOC_mesh.ID_jet_boundary],MOC_mesh.y[MOC_mesh.ID_jet_boundary],'go-')
# MOC_mesh.compute_thrust('nearest',10)
#plt.plot(MOC_mesh.x[MOC_mesh.ID_contour_chr],MOC_mesh.y[MOC_mesh.ID_contour_chr],'bo-')

#plt.plot(MOC_mesh.x[MOC_mesh.ID_contour_chr[-1]],0,'bo')
#plt.plot(spike.x,interpolate.splev(spike.x,tck,der=0),spike.lip_x,spike.lip_y,'X')
#plt.axis('equal')
# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.scatter(MOC_mesh.x,MOC_mesh.y,c=MOC_mesh.M, cmap = cm.coolwarm)
# #plt.plot(MOC_mesh.x,MOC_mesh.y,'r.')
# # cax = ax.imshow(MOC_mesh,interpolation='nearest',cmap=cm.coolwarm)

# # cbar = fig.colorbar(cax)
# plt.colorbar(ax=ax)
#plt.axis('equal')
plt.show()
