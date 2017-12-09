from plug_nozzle_angelino import *
import matplotlib.pyplot as plt
from scipy import interpolate
from scipy import optimize
import numpy as np



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


spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000)

tck = interpolate.splrep(spike.x,spike.y,full_output=0)

line = lambda x: np.tan(-np.pi/4)*(x-spike.lip_x) + spike.lip_y


## METHOD FOR CALCULATING FUNCTION INTERCEPTS
x = optimize.brentq(lambda x: interpolate.splev(x,tck,der=0) - line(x),spike.x[0],spike.x[-1]) 
 

plt.plot(spike.x,spike.y,spike.lip_x,spike.lip_y,'X',x,line(x),'o')
plt.plot(np.array([spike.lip_x,x]),line(np.array([spike.lip_x,x])))

#plt.plot(spike.x,interpolate.splev(spike.x,tck,der=1))
plt.show()

plt.axis('equal')
#### MOC DESCRIPTION

class chr_mesh():
    def __init__(self):
        pass

    # important variables: x, y, theta, mu, M, T, V, V_l, rho


    # two lists of point indexes:
    #   1 for left running characteristics, 1 for right running characteristics

    ### ALGORITHM STEPS

    # 1. Create initial data line:
    #   1.1 Calculate maximum expansion ratio for altitude
    #   1.2 Initialize Mach array for initial data line
