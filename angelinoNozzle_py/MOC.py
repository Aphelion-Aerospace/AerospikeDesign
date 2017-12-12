from plug_nozzle_angelino import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import optimize
import numpy as np
import gasdynamics as gd 



#### NOZZLE INITIAL CONDITIONS
n = 20

r_e = 0.067/2 #0.034 # likely too large
expansion_ratio = 6.64 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2343# np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 

p_inf = 1 # bar


#### DESIGN OF SPIKE

spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000)

tck = interpolate.splrep(spike.x,spike.y,full_output=0)

line = lambda x: np.tan(-np.pi/3)*(x-spike.lip_x) + spike.lip_y


## METHOD FOR CALCULATING FUNCTION INTERCEPTS
x = optimize.brentq(lambda x: interpolate.splev(x,tck,der=0) - line(x),spike.x[0],spike.x[-1]) 
 

plt.plot(spike.x,spike.y,spike.lip_x,spike.lip_y,'X',x,line(x),'o')
#plt.plot(np.array([spike.lip_x,x]),line(np.array([spike.lip_x,x])))
plt.axis('equal')
#plt.plot(spike.x,interpolate.splev(spike.x,tck,der=1))
#plt.show()


#### MOC DESCRIPTION
slope_init = -(spike.lip_x-spike.x[0])/(spike.lip_y-spike.y[0]); 

class chr_point():
    def __init__(self,gamma,x,y,theta,W):
        self.x = x; self.y = y; self.theta = theta; self.W = W
        self.mu = np.arcsin(np.arctan((gamma-1)/2*(1/self.W**2-1)))
        self.M = 1/np.sin(self.mu)
        # will calculate other properties later (p etc.) to vectorize
class chr_mesh():
    def __init__(self,spike,gamma,altitude,n):
        # constants of iteration
        tck = interpolate.splrep(spike.x,spike.y,full_output=0)
        (p_atm,T_atm,rho_atm) = gd.standard_atmosphere(altitude)
        self.gamma = gamma
        self.PR = spike.p_c/(p_atm*10**-5)

        # Point storage
        self.chr_dict = {} # dictionary for storing points based on ID
        self.chr_array = np.array([]) # array for storing points based on ID

        self.ID_not_computed =[] 
        # construction of initial expansion fan
        M_max = gd.PR_expansion_mach(self.PR,gamma)

        mach_fan = np.linspace(1.01,M_max,n)

        theta_fan = gd.prandtl_meyer(mach_fan,gamma)

        angle_fan = gd.mach_angle(mach_fan,degrees=1)



        # TODO: COMPUTE EXPANSION FAN POINTS (initial data line)
        # TODO: COMPUTE CHAR. NET DOWN TO NOZZLE BASE (when A.x > spike.x[-1], switch centre point jet boundary)
        ## after computing initial fan, load all point id's into ordered list showing right and left running characteristics that have not been used to compute a downstream point
        # TODO: COMPUTE NOZZLE BASE PRESSURE
        # TODO: COMPUTE PLUME DOWNSTREAM OF BASE
        # TODO: COMPUTE THRUST PRODUCED
        print(angle_fan)


    def general_point(A,B):
        # given points A & B (B.y>A.y) in char mesh, computes 3rd point
        x_c = (A.x*np.tan(A.theta+A.mu)-A.y+B.y-B.x*np.tan(B.theta-B.mu))/(np.tan(A.theta+A.mu)-np.tan(B.theta-B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta+A.mu) + A.y
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta-B.mu)
        theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A/A.y*(x_c-A.x))+B.W+B.W*(B.theta*np.tan(B.mu)+m_B/B.y*(x_c-B.x)))/(A.W*np.tan(A.mu)+B.W*np.tan(B.mu))
        W_c = A.W + A.W(np.tan(A.mu(theta_c-A.theta))+l_A/A.y*(x_c-A.x))
        self.chr_dict[self.max_ID + 1] = chr_point(self.gamma,x,y,theta,W)
    def same_fam_point(A,B):
        x_c = (A.x*np.tan(A.theta+A.mu)-A.y+B.y-B.x*np.tan(B.theta-B.mu))/(np.tan(A.theta+A.mu)-np.tan(B.theta-B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta+A.mu) + A.y
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta+B.mu)
        theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A/A.y*(x_c-A.x))+B.W+B.W*(B.theta*np.tan(B.mu)+m_B/B.y*(x_c-B.x)))/(A.W*np.tan(A.mu)+B.W*np.tan(B.mu))
        W_c = A.W + A.W(np.tan(A.mu(theta_c-A.theta))+l_A/A.y*(x_c-A.x))
        self.chr_dict[self.max_ID + 1] = chr_point(self.gamma,x,y,theta,W)

    def centre_point(A,B):
        pass

    def jet_boundary_point(A,B):
        pass
    # important variables: x, y, theta, mu, M, T, V, V_l, rho

    # two lists of point indexes:
    #   1 for left running characteristics, 1 for right running characteristics

    ### ALGORITHM STEPS



chr_mesh(spike,gamma,10000,10)

# tan_x = np.array([spike.x[0],spike.x[-1]])
# tan_y = slope_init*(tan_x-spike.lip_x)+spike.lip_y

# tan_yi = -1/slope_init*(tan_x-spike.lip_x)+spike.lip_y

# plt.plot(tan_x,tan_y,tan_x,tan_yi)


M = optimize.brentq(lambda M: gd.expansion_ratio_zero(1,M,gamma,expansion_ratio),1,10)

# print(M)
# print(spike.M[-1])
# print(gd.prandtl_meyer(spike.M[-1],gamma,degrees=1))
# print(np.arctan(slope_init)*180/np.pi)
#plt.show()
