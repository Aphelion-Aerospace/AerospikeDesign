import numpy as np 
from scipy import optimize
import matplotlib.pyplot as plt
from matplotlib import cm 
from mpl_toolkits.mplot3d import Axes3D
import gasdynamics as gd

def calc_line_intr(m1,x1,y1,m2,x2,y2):
    # takes in parameters for two lines in point slope form and returns their intercept
    if abs(m1) > 10**12:
        m1 = np.inf
    if abs(m2) > 10**12:
        m2 = np.inf
     
    if(np.isinf(m1)):
        if(np.isinf(m2)):
            raise ValueError('Both slopes inf, no intercept')
        else:
            x = x1
            y = m2*(x-x2) + y2
    elif(np.isinf(m2)):
        x = x2
        y = m1*(x-x1) + y1
    else:
        x = (y2 - y1 - m2*x2 + m1*x1)/(m1-m2)
        y = m1*(x-x1) + y1
    
    return x,y

class bell_nozzle:
	def __init__(self,expansion_ratio,A_t,gamma,T_c,p_c,n):
	
		self.A_t = A_t
		self.gamma = gamma
		self.n = n

		self.r_t = np.sqrt(A_t/np.pi)
		#self.r_t = 0.07
		self.no_points = int(1/2*n*(n+3))
		self.M_e = optimize.fsolve(lambda M: gd.expansion_ratio_zero(1,M,gamma,expansion_ratio),5) # should be edited to relate to expansion ratio

		self.K_n = np.zeros(self.no_points)

		self.K_p = np.zeros(self.no_points)

		self.theta = np.zeros(self.no_points)

		self.nu = np.zeros(self.no_points)

		self.M = np.zeros(self.no_points)

		self.mu = np.zeros(self.no_points)

		self.X = np.zeros(self.no_points)

		self.Y = np.zeros(self.no_points)

		self.theta_max = gd.prandtl_meyer(self.M_e,self.gamma)/2

		self.init_data_line()

		self.bell_fan_completion()

		T_ratio,p_ratio,rho_ratio,a_ratio = gd.isentropic_ratios(0,self.M,self.gamma)
		
		self.T = T_c*T_ratio
		self.P = p_c*p_ratio

		self.truncate_nozzle()
	def bell_fan_completion(self):
		exp_fan = self.n-1
		count = 0
		for i in range(self.n+1,self.no_points):
			if (count == 0):
				# centre line cross
				#print(str(i) + ' centre line')
				count += 1
				self.theta[i] = 0
				self.K_n[i] = self.K_n[i-exp_fan-1]
				self.nu[i] = self.K_n[i]
				self.K_p[i] = -self.nu[i]
				self.M[i] = optimize.fsolve(lambda M: gd.prandtl_meyer_zero(M,self.nu[i],self.gamma),1.1)
				self.mu[i] = gd.mach_angle(self.M[i])
				self.X[i],self.Y[i] = calc_line_intr(np.tan(self.theta[i-exp_fan-1]-self.mu[i-exp_fan-1]),self.X[i-exp_fan-1],self.Y[i-exp_fan-1],0,0,0) 

			elif (count == exp_fan):
				# side
				count = 0
				exp_fan -= 1
				self.theta[i] = self.theta[i-1]
				self.nu[i] = self.nu[i-1]
				self.M[i] = self.M[i-1]
				self.mu[i] = self.mu[i-1]
				self.K_n[i] = self.K_n[i-1]
				self.K_p[i] = self.K_p[i-1]				
				self.X[i],self.Y[i] = calc_line_intr(np.tan(self.theta[i-1] + self.mu[i-1]),self.X[i-1],self.Y[i-1],np.tan(1/2*(self.theta[i]+self.theta[i-exp_fan-2])),self.X[i-exp_fan-2],self.Y[i-exp_fan-2])
			else:
				# mid points
				#print(str(i) + ' mid points')
				count += 1
				self.K_n[i] = self.K_n[i-exp_fan-1]
				self.K_p[i] = self.K_p[i-1]
				self.theta[i] = 1/2*(self.K_n[i]+self.K_p[i])
				self.nu[i] = 1/2*(self.K_n[i]-self.K_p[i])
				self.M[i] = optimize.fsolve(lambda M: gd.prandtl_meyer_zero(M,self.nu[i],self.gamma),1.1)
				self.mu[i] = gd.mach_angle(self.M[i])
				self.X[i],self.Y[i] = calc_line_intr(np.tan(self.theta[i-exp_fan-1]-self.mu[i-exp_fan-1]),self.X[i-exp_fan-1],self.Y[i-exp_fan-1],np.tan(self.theta[i-1] + self.mu[i-1]),self.X[i-1],self.Y[i-1])		

	def init_data_line(self):
		# must be edited to remove hard coding
		
		self.theta[0:self.n] = np.linspace(0.0001,self.theta_max,self.n)
		#self.theta = self.theta*np.pi/180

		self.nu[0:self.n] = self.theta[0:n]

		self.M[0:self.n] = optimize.fsolve(lambda M: gd.prandtl_meyer_zero(M,self.nu[0:self.n],self.gamma),1.1*np.ones(self.n))

		self.mu[0:self.n] = gd.mach_angle(self.M[0:self.n])

		self.X[0],self.Y[0] = calc_line_intr(np.tan(self.theta[0] - self.mu[0]),0,self.r_t,0,0,0)

		for i in range(1,n):
			self.X[i],self.Y[i] = calc_line_intr(np.tan(self.theta[i] - self.mu[i]),0,self.r_t,np.tan(self.theta[i-1] + self.mu[i-1]),self.X[i-1],self.Y[i-1])

		self.theta[n] = self.theta[n-1]
		self.nu[n] = self.nu[n-1]
		self.M[n] = self.M[n-1]
		self.mu[n] = self.mu[n-1]

		self.K_n = self.theta+self.nu
		self.K_p = self.theta-self.nu

		self.X[n],self.Y[n] = calc_line_intr(np.tan(1/2*(self.theta_max + self.theta[n])),0,self.r_t,np.tan(self.theta[n-1]+self.mu[n-1]),self.X[n-1],self.Y[n-1])

	def truncate_nozzle(self):
		idx = self.X <= self.X.max()*0.7
		self.Y = self.Y[idx]; self.X = self.X[idx]; self.M = self.M[idx]; self.theta = self.theta[idx]; self.K_n = self.K_n[idx]; self.K_p = self.K_p[idx]; self.mu = self.mu[idx]; self.T = self.T[idx]; self.P = self.P[idx]

# END OF HELPER FUNCTIONS #############################

expansion_ratio = 8.1273
A_t = 0.00014
p_c = 34.474
T_c = 2833.63
gamma = np.mean([1.2534,1.2844])

n = 150 #no. of expansion waves

no_points = 1/2*n*(n+3)


b1 = bell_nozzle(expansion_ratio,A_t,gamma,T_c,p_c,n)

#plt.plot(0,0,'rx',0,b1.r_t,'rx',b1.X,b1.Y,'.',markersize=1)
#plt.contour(b1.X,b1.Y,b1.P,200)
#print(b1.T)

#print(b1.T[0])

plt.scatter(b1.X,b1.Y,c=b1.T,cmap=cm.coolwarm)
plt.colorbar()
print(b1.X.max())
plt.axis('equal')
plt.show()