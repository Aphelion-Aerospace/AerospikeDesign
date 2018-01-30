import numpy as np 
import gasdynamics as gd 
import matplotlib.pyplot as plt
from scipy import optimize
from scipy import interpolate
from matplotlib import cm 
import os

class plug_nozzle:
	def __init__(self,expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 1):
		# input design parameters	
		self.expansion_ratio = expansion_ratio
		self.A_t = r_e**2*np.pi/expansion_ratio
		self.r_e = r_e
		self.gamma = gamma
		self.n = n
		self.truncate_ratio = truncate_ratio

		self.T_c = T_c
		self.p_c = p_c
		self.a_c = a_c
		self.rho_c = rho_c

		# calculated design parameters
		self.A_e = self.A_t*self.expansion_ratio
		self.r_b = np.sqrt(-self.A_e/np.pi + self.r_e**2)
		self.M_e = optimize.fsolve(lambda M: gd.expansion_ratio_zero(1,M,self.gamma,self.expansion_ratio),5)

		# DESIGN OF NOZZLE, FUNCTION ORDER IS IMPORTANT
		# NON-OPTIONAL FUNCTION RUNS
		self.design_nozzle()

		self.truncate_nozzle()

		self.calc_flow_properties()

		self.arc_length_coord()

		# OPTIONAL FUNCTION CONSTANTS
		self.converge_section = 0 # whether the converging section has been designed


	## NON-OPTIONAL FUNCTIONS
	def design_nozzle(self):	
		# discrete contour design variables
		self.M = np.linspace(1,self.M_e,self.n)	
		self.A = self.A_t*gd.expansion_ratio(1,self.M,self.gamma)
		self.alpha = gd.prandtl_meyer(self.M_e,self.gamma) - gd.prandtl_meyer(self.M,self.gamma) + gd.mach_angle(self.M)
		self.l = (self.r_e - np.sqrt(np.abs(self.r_e**2 - (self.A*self.M*np.sin(self.alpha)/np.pi))))/np.sin(self.alpha)

		self.x = self.l*np.cos(self.alpha)
		self.y = self.l*np.sin(self.alpha)

		self.centre_spike()
		
		self.length = self.x.max()

	def centre_spike(self):
		self.lip_x = -self.x.min()
		self.lip_y = self.r_e 

		self.x = self.x - self.x.min()
		self.y = self.r_e - self.y


	def truncate_nozzle(self):
		# based on: Marcello Onofri, "Plug Nozzles: Summary of Flow Features and Engine Performance", University of Rome, 01 Jan 2006, American Institue of Aeronautics and Astronautics

		# Truncating to about 20% of the original length will produce and efficiency of of 0.82-0.97 for a pressure ratio of 8.9-200 (as opposed to 0.98-0.99 for full length nozzle)
		idx = self.x <= self.x.max()*self.truncate_ratio#x.max()#*0.2
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

	## OPTIONAL-FUNCTIONS
	def update_contour(self,x,y,M,x_centre_spike=0):
		# Updates the spike contour with new (x,y) points with known values of M at each point
		self.x = x; self.y = y; self.M = M

		# optionally centre spike about x-axis
		if(x_centre_spike):
			self.centre_spike()

		# update flow properties based on isentropic expansion
		self.calc_flow_properties()

		# update arc length coordinates
		self.arc_length_coord()

		# update exit mach number
		self.M_e = M[-1]

		# update expansion ratio
		self.expansion_ratio = gd.expansion_ratio(1,self.M_e)

		self.A_t = self.r_e**2*np.pi/self.expansion_ratio
		print("Warning, throat area update not complete, assumes perfect isentropic expansion from throat to exit")
		
		# update area estimation 
		self.A = self.A_t*gd.expansion_ratio(1,self.M,self.gamma)
		# update base radius
		self.r_b = self.y[-1]

		# update exit area
		self.A_e = np.pi*(self.r_e**2-self.r_b**2)

		if(self.converge_section):
			print("Warning, congerence section not updated. Run self.converge_section(args) again to define a new convergence section.")


	def calc_ideal_thrust(self,p_atm):
		# calculates ideal thrust
		# F = m_dot*V_e + (P_e - P_o)A_e
		p_e = self.p[-1]

		thrust = self.rho[0]*self.V[0]*self.A_t*self.V[-1] + (p_e-p_atm)*self.A_e
		return thrust 

	def define_compression(self,r1,r2,slope,conv_length,n):
		self.converge_section = 1
		tck = interpolate.splrep(self.x,self.y)

		alpha = np.arctan(-1/interpolate.splev(self.x[0],tck,der=1))

		x1 = -r1*np.cos(alpha); x2 = x1 

		y1 = self.y[0]-r1*np.sin(alpha)

		y2 = r1 + y1 - r2

		beta = np.arctan(-1/slope)+np.pi

		x_str_bnd = x2 + r2*np.cos(beta)

		y_str_bnd = y2 + r2*np.sin(beta)

		def conv_geom(x):
			if (x > x1):
				theta = np.arccos((x-x1)/r1)
				y = r1*np.sin(theta) + y1
			elif (x > x_str_bnd):
				theta = np.arccos((x-x2)/r2)
				y = r2*np.sin(theta) + y2
			else:
				y = slope*(x-x_str_bnd) + y_str_bnd

			return y

		x_init = x_str_bnd - np.sqrt(conv_length**2/(1+slope**2))

		self.conv_x = np.linspace(x_init,self.x[0],n)
		self.conv_y = np.ones(self.conv_x.shape)
		for i in range(len(self.conv_x)):
			self.conv_y[i] = conv_geom(self.conv_x[i])
		#print(self.conv_x)

	def plot_contour(self,ax):

		if (self.converge_section):
			ax.plot(self.conv_x,self.conv_y)

		ax.plot(self.x,self.y)
		ax.plot(self.lip_x,self.lip_y,'rx')
		ax.plot(self.x,np.zeros(self.x.shape),'k--')

	def save_to_csv(self):
		if not os.path.exists('plug_csv'):
			os.makedirs('plug_csv')

		csv_array = np.array([self.x,self.y,self.s,self.p,self.T,self.M,self.A,self.a,self.V,self.rho])
		np.savetxt('plug_csv/aerospike_diverge_contour.csv', csv_array.T, delimiter = ',')
		with open('plug_csv/aerospike_diverge_contour.csv','r') as original: data = original.read()
		with open('plug_csv/aerospike_diverge_contour.csv','w') as modified: modified.write('x,y,s,p,T,M,A,a,V,rho\n' + data)

		csv_array = np.array([[self.lip_x],[self.lip_y]])		
		np.savetxt('plug_csv/aerospike_lip_coordinates.csv',csv_array.T,delimiter =',')
		with open('plug_csv/aerospike_lip_coordinates.csv','r') as original: data = original.read()
		with open('plug_csv/aerospike_lip_coordinates.csv','w') as modified: modified.write('lip x,lip y\n' + data)

		csv_array = np.array([self.conv_x,self.conv_y])
		np.savetxt('plug_csv/aerospike_converge_contour.csv', csv_array.T,delimiter = ',')
		with open('plug_csv/aerospike_converge_contour.csv','r') as original: data = original.read()
		with open('plug_csv/aerospike_converge_contour.csv','w') as modified: modified.write('Converging x,Converging y\n' + data)

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