import numpy as np 
import gasdynamics as gd 
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate
from matplotlib import cm 

class plug_nozzle:
	def __init__(self,expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n):
		# input design parameters	
		self.expansion_ratio = expansion_ratio
		self.A_t = A_t
		self.r_e = r_e
		self.gamma = gamma
		self.n = n

		self.T_c = T_c
		self.p_c = p_c
		self.a_c = a_c
		self.rho_c = rho_c

		# calculated design parameters
		self.A_e = self.A_t*self.expansion_ratio
		self.r_b = np.sqrt(-self.A_e/np.pi + self.r_e**2)
		self.M_e = optimize.fsolve(lambda M: gd.expansion_ratio_zero(1,M,self.gamma,self.expansion_ratio),5)

		# DESIGN OF NOZZLE, FUNCTION ORDER IS IMPORTANT

		self.design_nozzle()

		self.truncate_nozzle()

		self.calc_flow_properties()

		self.arc_length_coord()

	def design_nozzle(self):	
		# discrete contour design variables
		self.M = np.linspace(1,self.M_e,self.n)	
		self.A = self.A_t*gd.expansion_ratio(1,self.M,self.gamma)
		self.alpha = gd.prandtl_meyer(self.M_e,self.gamma) - gd.prandtl_meyer(self.M,self.gamma) + gd.mach_angle(self.M)
		self.l = (self.r_e - np.sqrt(np.abs(self.r_e**2 - (self.A*self.M*np.sin(self.alpha)/np.pi))))/np.sin(self.alpha)

		self.x = self.l*np.cos(self.alpha)
		self.y = self.l*np.sin(self.alpha)

		self.lip_x = -self.x.min()
		self.lip_y = self.r_e 

		self.x = self.x - self.x.min()
		self.y = self.r_e - self.y
		
		self.length = self.x.max()

	def truncate_nozzle(self):
		# based on: Marcello Onofri, "Plug Nozzles: Summary of Flow Features and Engine Performance", University of Rome, 01 Jan 2006, American Institue of Aeronautics and Astronautics

		# Truncating to about 20% of the original length will produce and efficiency of of 0.82-0.97 for a pressure ratio of 8.9-200 (as opposed to 0.98-0.99 for full length nozzle)
		idx = self.x <= self.x.max()*0.2
		self.M = self.M[idx]; self.A = self.A[idx]; self.alpha = self.alpha[idx]; self.l = self.l[idx]; self.x= self.x[idx]; self.y = self.y[idx]; 

	def calc_flow_properties(self):
		T_ratio,p_ratio,rho_ratio,a_ratio = gd.isentropic_ratios(0,self.M,self.gamma)
		self.T = self.T_c*T_ratio
		self.p = self.p_c*p_ratio	
		self.a = self.a_c*a_ratio
		self.V = self.a*self.M	
		self.rho = self.rho_c*rho_ratio

	def arc_length_coord(self):
		y_dummy = self.y[1:] - self.y[:-1]
		x_dummy = self.x[1:] - self.x[:-1]
		s_dummy = np.sqrt(y_dummy**2+x_dummy**2)

		s_dummy = np.concatenate((np.array([0]),s_dummy))

		self.s = np.zeros(s_dummy.shape)

		for i in range(1,s_dummy.shape[0]):
			self.s[i] = s_dummy[i] + self.s[i-1]

	def internal_compression(self):
		pass


###
# End of helper function / class descriptions
###

#design for 30,000

# r_e = 0.072/2 #0.034 # likely too large
# expansion_ratio = 6.64 #8.1273
# A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
# gamma = 1.2343# np.mean([1.2534,1.2852])
# T_c = 2833.63
# p_c = 34.474
# rho_c = 3.3826
# a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 

# print('Sound speed: ' + str(a_c))

# plug1 = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000)

# plt.plot(plug1.x,plug1.y, label='Aerospike Contour')#c=plug1.rho,cmap=cm.coolwarm)
# plt.plot(plug1.lip_x,plug1.lip_y,'rx',label='Lip Location')
# #plt.colorbar()
# plt.plot([0,plug1.x.max()],[0,0], 'k--',label='Centre Line')
# plt.legend()
# print('Distance above r_t: ' + str(plug1.lip_y - plug1.y[0]))
# plt.xlabel('x (m)')
# plt.ylabel('r (m)')

# m = (plug1.lip_y - plug1.y[0])/(plug1.lip_x - plug1.x[0])

# m = -1/m

# print('Flow angle at throat: ' + str(180/np.pi*np.tan(m)-180))

# max_y = m*(-plug1.lip_x) + plug1.lip_y

# # plt.plot(0,max_y,'gx')
# plt.axis('equal')



# print('radius of curvature near the throat: ' + str(2*np.sqrt((plug1.lip_x - plug1.x[0])**2 + (plug1.lip_y - plug1.y[0])**2)))

# csv_array = np.array([plug1.x,plug1.y,plug1.s,plug1.p,plug1.T,plug1.M,plug1.A,plug1.a,plug1.V,plug1.rho])

# np.savetxt('aerospike_contour.csv', csv_array.T, delimiter = ',')

# ## plots of p,T,M,a,V,rho

# fig1, ((ax1,ax2,ax3),(ax4,ax5,ax6)) = plt.subplots(2,3)

# ax1.plot(plug1.x*100,plug1.p*100)
# #ax1.set_xlabel('x (cm)')
# ax1.set_ylabel('kPa')
# ax1.set_title('Pressure on Contour Surface')
# ax1.grid()

# ax2.plot(plug1.x*100,plug1.T)
# #ax2.set_xlabel('x (cm)')
# ax2.set_ylabel('K')
# ax2.set_title('Temperature on Contour Surface')
# ax2.grid()

# ax3.plot(plug1.x*100,plug1.M)
# #ax3.set_xlabel('x (cm)')
# ax3.set_ylabel('M')
# ax3.set_title('Mach on Contour Surface')
# ax3.grid()

# ax4.plot(plug1.x*100,plug1.a)
# ax4.set_xlabel('x (cm)')
# ax4.set_ylabel('m/s')
# ax4.set_title('Sound Speed on Contour Surface')
# ax4.grid()

# ax5.plot(plug1.x*100,plug1.V)
# ax5.set_xlabel('x (cm)')
# ax5.set_ylabel('m/s')
# ax5.set_title('Velocity on Contour Surface')
# ax5.grid()

# ax6.plot(plug1.x*100,plug1.rho)
# ax6.set_xlabel('x (cm)')
# ax6.set_ylabel('KG/CU')
# ax6.set_title('Density on Contour Surface')
# ax6.grid()

# plt.show()