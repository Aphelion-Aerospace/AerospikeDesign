from plug_nozzle_angelino import *
import matplotlib.pyplot as plt
from scipy import interpolate
import numpy as np

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

n = 20


r_e = 0.067/2 #0.034 # likely too large
expansion_ratio = 6.64 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2343# np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 


spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000,truncate_ratio = 0.2)

#tck = interpolate.splrep(spike.x,spike.y)



plt.plot(spike.x,spike.y,spike.lip_x,spike.lip_y,'X')

# defining converging section

m = (spike.lip_y  - spike.y[0])/(spike.lip_x - spike.x[0])

theta_init = np.arctan(m)

x_int,y_int = calc_line_intr(m,spike.x[0],spike.y[0],np.inf,-0.008,0)

pts = 100

r_inner = np.linspace(np.sqrt((spike.x[0]-x_int)**2 + (spike.y[0]-y_int)**2),np.sqrt((spike.x[0]-x_int)**2 + (spike.y[0]-y_int)**2),pts)

r_outer = np.linspace(np.sqrt((spike.lip_y-y_int)**2 + (spike.lip_x-x_int)**2),np.sqrt((spike.lip_y-y_int)**2 + (spike.lip_x-x_int)**2)+0.002,pts)

theta = np.linspace(theta_init,np.pi/2,pts)

x_outer = r_outer*np.cos(theta) + x_int
y_outer = r_outer*np.sin(theta) + y_int

x_inner = r_inner*np.cos(theta) + x_int
y_inner = r_inner*np.sin(theta) + y_int

csv_array = np.array([spike.x,spike.y,spike.s,spike.p,spike.T,spike.M,spike.A,spike.a,spike.V,spike.rho])

np.savetxt('aerospike_contour.csv', csv_array.T, delimiter = ',')

csv_array = np.array([x_outer,y_outer,x_inner,y_inner])

np.savetxt('converging_contour.csv',csv_array.T, delimiter = ',')



plt.plot(x_outer,y_outer,x_inner,y_inner)

plt.axis('equal')

plt.show()
