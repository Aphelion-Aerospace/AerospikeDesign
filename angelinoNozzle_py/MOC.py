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
gamma = 1.2381 #np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474
rho_c = 3.3826
R = (1-1/gamma)*1.8292*1000#1.8292 1.6196
a_c = np.sqrt(gamma*R*T_c)

p_inf = 1 # bar


#### DESIGN OF SPIKE

spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,10000,truncate_ratio =0.2)

spike.y = spike.y*-1
spike.lip_y = spike.lip_y*-1

line = lambda x: np.tan(-np.pi/3)*(x-spike.lip_x) + spike.lip_y


## METHOD FOR CALCULATING FUNCTION INTERCEPTS

tck = interpolate.splrep(spike.x,spike.y,full_output=0)
plt.plot(spike.x,interpolate.splev(spike.x,tck,der=0),spike.lip_x,spike.lip_y,'X')
#plt.plot(np.array([spike.lip_x,x]),line(np.array([spike.lip_x,x])))
plt.axis('equal')
#plt.show()
#plt.plot(spike.x,interpolate.splev(spike.x,tck,der=1))

#### MOC DESCRIPTION

class chr_point():
    def __init__(self,gamma,x,y,theta,W):
        self.x = x; self.y = y; self.theta = theta; self.W = W
        self.mu = np.arcsin(np.sqrt((gamma-1)/2*(1/W**2-1)))
        # print('mu! + ' + str(self.mu))
        self.M = 1/np.sin(self.mu)
        # will calculate other properties later (p etc.) to vectorize
class chr_mesh():
    def __init__(self,spike,gamma,altitude,n,downstream_distance,plot_chr=0):

        self.spike =spike; self.gamma = gamma; self.altitude =altitude; self.n = n; self.downstream_distance = downstream_distance
        # constants of iteration

        self.plot_chr = plot_chr

        self.slope_init = np.arctan(-(spike.lip_x-spike.x[0])/(spike.lip_y-spike.y[0])); 
        print(self.slope_init*180/np.pi)

        self.tck = interpolate.splrep(spike.x,spike.y,full_output=0)
        (self.p_atm,self.T_atm,self.rho_atm) = gd.standard_atmosphere(altitude)
        self.gamma = gamma
        self.PR = self.spike.p_c/(self.p_atm*10**-5)
        self.V_l = np.sqrt(2/(self.gamma-1))*spike.a_c

        # Point storage

        self.chr_array = np.array([]) # array for storing points based on ID

        # logs characteristics that have not been intercepted
        self.ID_left_chr = [] 
        self.ID_right_chr = []
        self.ID_jet_boundary = []
        self.ID_contour_chr = []
        #self.ID_next_chr_jet = []
        #self.ID_compression_chr = []

        self.ID = 0 # ID of next point, not yet computed 
        # COMPUTE EXPANSION FAN POINTS (initial data line)
        # construction of initial expansion fan
       
        self.compute_initial_expansion_fan()
        # TODO: COMPUTE CHAR. NET DOWN TO NOZZLE BASE (when A.x > spike.x[-1], switch centre point jet boundary)
        ## after computing initial fan, load all point id's into ordered list showing right and left running characteristics that have not been used to compute a downstream point
        # TODO: COMPUTE NOZZLE BASE PRESSURE
        # while (self.chr_array[self.ID_contour_chr[0]].x <= spike.x.max()):
        #     pass
        # TODO: COMPUTE PLUME DOWNSTREAM OF BASE
        # TODO: COMPUTE THRUST PRODUCED 
        # for point in self.chr_array:
        #     plt.plot(point.x,point.y,'rX')       
        print(self.ID_left_chr)
        print(self.ID_right_chr)
        print(self.ID_jet_boundary)
        print(self.ID_contour_chr)

        # base conditions
        self.p_b = self.base_pressure()




        self.contour_fan = 1
        self.new_fan = 1
        self.first_base_intercept = 1
        #while (self.chr_array[self.ID_contour_chr[-1]].x <= spike.x.max()):
        for i in range(700):
            ## CONTOUR FAN
            # if (self.contour_fan):
            ID_temp = self.ID_right_chr.pop(0)
            #print(ID_temp)
            if(self.new_fan):
                # intersection
                self.new_fan = 0
                if (self.on_nozzle_contour(self.chr_array[ID_temp])):
                    new_point = self.contour_point(self.chr_array[ID_temp],plot_chr=self.plot_chr)
                    self.chr_array = np.append(self.chr_array,new_point)
                    self.ID += 1
                    # first point
                    ID_temp = self.ID_right_chr.pop(0)
                    new_point = self.general_point(self.chr_array[ID_temp],self.chr_array[self.ID-1],plot_chr=self.plot_chr)
                    self.chr_array = np.append(self.chr_array,new_point)
                    self.ID += 1     
                else:
                    if (self.first_base_intercept):
                        self.first_base_intercept = 0 
                        first_base_point = self.chr_array[self.ID_contour_chr[-1]]
                        M_b = gd.PR_expansion_mach(self.spike.p_c/self.p_b,self.gamma)
                        theta_b = gd.prandtl_meyer(M_b) - gd.prandtl_meyer(first_base_point.M)
                        W_b = first_base_intercept.W

                        first_base_point = chr_point(self.gamma,self.spike.x[-1],self.spike.y[-1],theta_b,W_b)
                        
                        new_point = self.internal_jet_point(first_base_point,self.chr_array[ID_temp])
                        self.chr_array = np.append(self.chr_array,new_point)
                        self.ID += 1

                    else: 
                        new_point = self.internal_jet_point(self.chr_array[self.ID_contour_chr[-1]],self.chr_array[ID_temp])
                        self.chr_array = np.append(self.chr_array,new_point)
                        self.ID += 1
            elif(ID_temp==-1):
                # self.ID_next_chr_jet.append(self.ID_left_chr.pop(0))   
                # plt.plot(self.chr_array[self.ID_jet_boundary[-1]].x,self.chr_array[self.ID_jet_boundary[-1]].y,'gx')

                new_point = self.general_point(self.chr_array[self.ID_jet_boundary[-1]],self.chr_array[self.ID-1],plot_chr=self.plot_chr)

                self.chr_array = np.append(self.chr_array,new_point)
                self.ID_left_chr.append(self.ID)

                self.ID += 1
                new_point = self.jet_boundary_point(self.chr_array[self.ID_jet_boundary[-1]],self.chr_array[self.ID-1],plot_chr=self.plot_chr)

                self.chr_array = np.append(self.chr_array,new_point)
                #plt.plot(new_point.x,new_point.y,'gx')
                self.ID_jet_boundary.append(self.ID)
                self.ID +=1
                self.contour_fan = 0   
                self.add_break_ID()
                self.new_fan = 1      
            else:
                temp1 = self.same_fam_point(self.chr_array[self.ID_right_chr[0]],self.chr_array[ID_temp])
                temp2 = self.general_point(self.chr_array[ID_temp],self.chr_array[self.ID-1]) 
                if (temp1.x < temp2.x) and (temp1.x>self.chr_array[ID_temp].x) :
                    #self.ID_right_chr.pop(-1); self.ID_left_chr.pop(-1)
                    #self.compression_offset += 1
                    plt.plot(self.chr_array[self.ID_right_chr[0]].x,self.chr_array[self.ID_right_chr[0]].y,'rx',self.chr_array[ID_temp].x,self.chr_array[ID_temp].y,'rx')
                    plt.plot(temp1.x,temp1.y,'gx')

                    plt.plot([self.chr_array[ID_temp].x, temp1.x],[self.chr_array[ID_temp].y, temp1.y])
                    self.chr_array = np.append(self.chr_array,temp1)
                    
                    new_point = self.general_point(self.chr_array[-1],self.chr_array[self.ID-1],plot_chr=1) 
                    self.ID += 1
                else:
                    new_point = temp2
                    plt.plot([self.chr_array[ID_temp].x, temp2.x],[self.chr_array[ID_temp].y, temp2.y],'k',[self.chr_array[self.ID-1].x, temp2.x],[self.chr_array[self.ID-1].y, temp2.y],'k')                    
                self.ID += 1
                self.chr_array = np.append(self.chr_array,new_point)
                
            # ## JET BOUNDARY FAN
            # else:
            #     ID_temp = self.ID_right_chr.pop(0)
            #     if(self.new_fan):

            #         new_point = self.general_point(self.chr_array[ID_temp],self.chr_array[self.ID_contour_chr[-1]])
            #         self.chr_array = np.append(self.chr_array,new_point)
            #         self.ID += 1
            #         self.new_fan = 0
            #     elif (ID_temp == -1):

            #         new_point=self.jet_boundary_point(self.chr_array[self.ID_next_chr_jet[-1]],self.chr_array[ID_temp -1])
            #         self.chr_array = np.append(self.chr_array, new_point)
            #         self.ID += 1
            #         self.contour_fan = 1
            #         self.new_fan = 1
            #         self.add_break_ID()
            #     else:
            #         new_point = self.general_point(self.chr_array[ID_temp],self.chr_array[ID_temp-1]) 
            #         self.chr_array = np.append(self.chr_array,new_point)

            #         self.ID += 1                   


    def compute_initial_expansion_fan(self):
        M_max = gd.PR_expansion_mach(self.PR,self.gamma)
        # print('M_max: ' + str(M_max))
        mach_fan = np.linspace(1.1,M_max,self.n)

        (T_ratio,p_ratio,rho_ratio,a_ratio) = gd.isentropic_ratios(0,mach_fan,self.gamma)

        V_fan = a_ratio*spike.a_c*mach_fan

        W_fan = V_fan/self.V_l

        theta_fan = -gd.prandtl_meyer(mach_fan,gamma) + self.slope_init

        angle_fan = gd.mach_angle(mach_fan)

        # print(180/np.pi*np.arcsin(np.sqrt((gamma-1)/2*(1/W_fan**2-1))))
        # print(W_fan)

        # print(angle_fan*180/np.pi)

        x_fan = np.ones(angle_fan.shape)*spike.lip_x

        y_fan = np.ones(angle_fan.shape)*spike.lip_y

        #print(theta_fan*180/np.pi)
        # print(gd.mach_angle_velocity_ratio(gd.prandtl_meyer(2.3,gamma),0.3,gamma))       

        initial_point = self.contour_point(chr_point(self.gamma,x_fan[0],y_fan[0],theta_fan[0],W_fan[0]),plot_chr=self.plot_chr)
        self.ID += 1
        self.ID_contour_chr.pop(0)
        self.chr_array = np.append(self.chr_array,initial_point)

        for point in x_fan[1:-1]:
            temp_point = chr_point(self.gamma,x_fan[self.ID],y_fan[self.ID],theta_fan[self.ID],W_fan[self.ID])
            new_point = self.general_point(temp_point,self.chr_array[self.ID-1],plot_chr=self.plot_chr)
            # adding to arrays
            self.chr_array = np.append(self.chr_array,new_point)
            self.ID += 1

        first_jet = chr_point(self.gamma,x_fan[-1],y_fan[-1],theta_fan[-1],W_fan[-1])
        second_jet = self.jet_boundary_point(first_jet,self.chr_array[self.ID-1],plot_chr=self.plot_chr)
        self.chr_array = np.append(self.chr_array,second_jet)
        self.ID_jet_boundary.append(self.ID)
        self.ID += 1
        self.add_break_ID()

    def general_point(self,A,B,plot_chr=0):
        # given points A & B (B.y>A.y) in char mesh, computes 3rd point
        x_c = (A.x*np.tan(A.theta+A.mu)-A.y+B.y-B.x*np.tan(B.theta-B.mu))/(np.tan(A.theta+A.mu)-np.tan(B.theta-B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta+A.mu) + A.y
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta-B.mu)
        theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A/A.y*(x_c-A.x))+B.W+B.W*(B.theta*np.tan(B.mu)+m_B/B.y*(x_c-B.x)))/(A.W*np.tan(A.mu)+B.W*np.tan(B.mu))
        W_c = A.W + A.W*(np.tan(A.mu*(theta_c-A.theta))+l_A/A.y*(x_c-A.x))

        self.ID_left_chr.append(self.ID)
        self.ID_right_chr.append(self.ID)
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'k',[B.x,x_c],[B.y,y_c],'k')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c)

    def same_fam_point(self,A,B,plot_chr=0):
        # DEPRECATED, SHOULD BE CHECKED
        # given points A & B (B.y>A.y) in char mesh, computes 3rd point in the case that the interception is of the same family of points
        x_c = (A.x*np.tan(A.theta+A.mu)-A.y+B.y-B.x*np.tan(B.theta+B.mu))/(np.tan(A.theta+A.mu)-np.tan(B.theta+B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta+A.mu) + A.y
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta+B.mu)
        theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A/A.y*(x_c-A.x))+B.W+B.W*(B.theta*np.tan(B.mu)+m_B/B.y*(x_c-B.x)))/(A.W*np.tan(A.mu)+B.W*np.tan(B.mu))
        W_c = A.W + A.W*(np.tan(A.mu*(theta_c-A.theta))+l_A/A.y*(x_c-A.x))
        # if (self.plot_chr):
        #     plt.plot([A.x,x_c],[A.y,y_c],'r',[B.x,x_c],[B.y,y_c],'r')
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'k',[B.x,x_c],[B.y,y_c],'k')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c)

    def contour_point(self,A,plot_chr=0):
        # Given a chr_point A, computes a chr_point that intersects with the nozzle boundary
        line = lambda x: np.tan(A.theta+A.mu)*(x-A.x) + A.y
        # print('line: ' + str((A.theta+A.mu)*180/np.pi))
        intercept = lambda x: interpolate.splev(x,self.tck,der=0) - line(x)
        x_c = optimize.brentq(intercept,self.spike.x[0],spike.x[-1])
        y_c = line(x_c)  
        theta_c = np.arctan(interpolate.splev(x_c,self.tck,der=1)) 
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)  
        W_c = A.W + A.W*((theta_c-A.theta)*np.tan(A.mu)+l_A/A.y*(x_c-A.x))
        #plt.plot(x_c,y_c,'rX')
        self.ID_contour_chr.append(self.ID)
        if(plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'k')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c)

    def jet_boundary_point(self,A,B,plot_chr=0):
        # given points A and B (where A is the previous jet boundary point and A.x < B.x) return the next boundary point
        x_c = (A.x*np.tan(A.theta)-A.y+B.y-B.x*np.tan(B.theta-B.mu))/(np.tan(A.theta)-np.tan(B.theta-B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta) + A.y
        mu_c = A.mu
        W_c = A.W
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta-B.mu)
        theta_c = B.theta + (-(W_c-B.W)/B.W + m_B*(x_c-B.x)/B.y)/np.tan(B.mu)
        self.ID_jet_boundary.append(self.ID)
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'r',[B.x,x_c],[B.y,y_c],'r')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c)
    def internal_jet_point(self,A,B,plot_chr=0):
        # given points A and B (where A is previous internal jet boundary point)
        x_c = (A.x*np.tan(A.theta)-A.y+B.y - B.x*np.tan(B.theta + B.mu))/(np.tan(A.theta) - np.tan(B.theta + B.mu))
        y_c = (x_c - A.x)*np.tan(A.theta) + A.y
        W_c = A.W
        mu_c = A.mu 

        theta_c = B.theta + (-(W_c -B.W)/B.W + m_B*(x_c - B.x)/B.y)/np.tan(B.mu)
        self.ID_contour_chr.append(self.ID)
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'r',[B.x,x_c],[B.y,y_c],'r')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c)                

    def add_break_ID(self):
        #self.ID_left_chr.append(-1)
        self.ID_right_chr.append(-1)
        #self.ID_jet_boundary.append(-1)
        #self.ID_contour_chr.append(-1)      

    def base_pressure(self):
        # F = C_f *A_t * P_c => F = rho_t*A_t*V_t*V_e + (P_e-P_a)*A_e

        epsilon_b = np.pi*self.spike.x[-1]**2/self.spike.A_t

        M_eb = optimize.brentq(lambda M: gd.expansion_ratio_zero(1,M,self.gamma,(self.spike.A[-1]-np.pi*self.spike.y[-1]**2)/self.spike.A_t),1,self.spike.M[-1])

        (T_ratio,p_ratio,rho_ratio,a_ratio) = gd.isentropic_ratios(1,M_eb,self.gamma)

        V_eb = a_ratio*self.spike.a_c*M_eb

        p_eb = p_ratio*self.spike.p_c

        F_no_core = self.spike.rho[0]*self.spike.A_t*self.spike.V[0]*V_eb + (p_eb - self.p_atm)*(self.spike.A[-1]-np.pi*self.spike.y[-1]**2)

        C_f = F_no_core/(self.spike.A_t*self.spike.p_c)

        # print('Base pressure = ' + str(0.648*C_f/epsilon_b), ' Atm pressure = ' + str(self.p_atm))

        return self.p_atm

    def on_nozzle_contour(self,chr_point_obj):
        line = lambda x: np.tan(chr_point_obj.theta+chr_point_obj.mu)*(x-chr_point_obj.x) + chr_point_obj.y
        intercept = lambda x: interpolate.splev(x,self.tck,der=0) - line(x)
        return np.sign(intercept(self.spike.x[0])*intercept(self.spike.x[-1])) <= 0

    # important variables: x, y, theta, mu, M, T, V, V_l, rho

    # two lists of point indexes:
    #   1 for left running characteristics, 1 for right running characteristics

    ### ALGORITHM STEPS



MOC_mesh = chr_mesh(spike,gamma,100,30,spike.x.max(),plot_chr=1)

# tan_x = np.array([spike.x[0],spike.x[-1]])
# tan_y = slope_init*(tan_x-spike.lip_x)+spike.lip_y

# tan_yi = -1/slope_init*(tan_x-spike.lip_x)+spike.lip_y

# plt.plot(tan_x,tan_y,tan_x,tan_yi)

plt.axis('equal')
plt.show()


M = optimize.brentq(lambda M: gd.expansion_ratio_zero(1,M,gamma,expansion_ratio),1,10)

# print(M)
# print(spike.M[-1])
# print(gd.prandtl_meyer(spike.M[-1],gamma,degrees=1))
# print(np.arctan(slope_init)*180/np.pi)
#plt.show()
