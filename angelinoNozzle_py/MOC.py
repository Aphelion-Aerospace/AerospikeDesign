from plug_nozzle_angelino import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from scipy import interpolate
from scipy import optimize
import numpy as np
import gasdynamics as gd 
from matplotlib import cm 
import copy

#### MOC DESCRIPTION

class chr_point():
    def __init__(self,gamma,x,y,theta,W,pt_type):
        self.x = x; self.y = y; self.theta = theta; self.W = W;self.pt_type = pt_type
        with np.errstate(invalid = 'ignore'):
            #pt_type = general, jet, contour, or N/A (temp points)
            self.mu = np.arcsin(np.sqrt((gamma-1)/2*(1/W**2-1)))
            # print('mu! + ' + str(self.mu))
            self.M = 1/np.sin(self.mu)

    def print_point(self):
        print('x = ' + str(self.x) + ' y = ' + str(self.y) + ' theta = ' + str(self.theta) + ' W = ' + str(self.W) + ' mu = ' + str(self.mu) + ' M = ' +str(self.M))

    def plot_point(self):
        plt.plot(self.x,self.y,'bx')
        # will calculate other properties later (p etc.) to vectorize
class chr_mesh():
    def __init__(self,spike,gamma,altitude,n,downstream_factor=1.1,plot_chr=0,clean_mesh=1):

        self.spike =copy.copy(spike); self.gamma = gamma; self.altitude =altitude; self.n = n

        self.downstream_factor = downstream_factor # percentage down after mesh cross with centre to continue meshing
        # constants of iteration

        self.plot_chr = plot_chr


        self.flip_plug() # flips sign of plug if required

        self.slope_init = np.arctan(-(self.spike.lip_x-self.spike.x[0])/(self.spike.lip_y-self.spike.y[0])); 
        #print(self.slope_init*180/np.pi)

        self.tck = interpolate.splrep(self.spike.x,self.spike.y,full_output=0)
        (self.p_atm,self.T_atm,self.rho_atm) = gd.standard_atmosphere([altitude])
        self.gamma = gamma
        self.PR = self.spike.p_c/(self.p_atm)
        self.V_l = np.sqrt(2/(self.gamma-1))*self.spike.a_c

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
        # print(self.ID_left_chr)
        # print(self.ID_right_chr)
        # print(self.ID_jet_boundary)
        # print(self.ID_contour_chr)

        # base conditions
        self.p_b = self.base_pressure()




        self.contour_fan = 1
        self.new_fan = 1
        self.first_base_intercept = 1
        self.centre_line_intercept = 0
        self.END_SIM = 0
        #while (self.chr_array[self.ID_contour_chr[-1]].x <= spike.x.max()):
        while self.chr_point_less_zero() and self.contour_converge():
        ## TODO: COMPUTE EXPANSION FAN UNTIL POINT IS > spike.length in which case, remove from all tracking lists and do not add to chr_array
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
                        #print(self.spike.p_c/self.p_b)
                        M_b = gd.PR_expansion_mach(self.spike.p_c/self.p_b,self.gamma)
                        #print(M_b)
                        theta_b = gd.prandtl_meyer(M_b,self.gamma) - gd.prandtl_meyer(first_base_point.M,self.gamma)
                        W_b = first_base_point.W

                        first_base_point = chr_point(self.gamma,self.spike.x[-1],self.spike.y[-1],theta_b,W_b,'contour')
                        new_point = self.internal_jet_point(first_base_point,self.chr_array[ID_temp],plot_chr=self.plot_chr)
                        self.chr_array = np.append(self.chr_array,new_point)
                        self.ID += 1

                    else: 
                        new_point = self.internal_jet_point(self.chr_array[self.ID_contour_chr[-1]],self.chr_array[ID_temp],plot_chr=self.plot_chr)
                        self.chr_array = np.append(self.chr_array,new_point)

                        self.ID += 1
            elif(ID_temp==-1):
                # self.ID_next_chr_jet.append(self.ID_left_chr.pop(0))   
                # plt.plot(self.chr_array[self.ID_jet_boundary[-1]].x,self.chr_array[self.ID_jet_boundary[-1]].y,'gx')

                new_point = self.general_point(self.chr_array[self.ID_jet_boundary[-1]],self.chr_array[self.ID-1],plot_chr=self.plot_chr)

                self.chr_array = np.append(self.chr_array,new_point)
                #self.ID_left_chr.append(self.ID)

                self.ID += 1
                new_point = self.jet_boundary_point(self.chr_array[self.ID_jet_boundary[-1]],self.chr_array[self.ID-1],plot_chr=self.plot_chr)

                self.chr_array = np.append(self.chr_array,new_point)
                #plt.plot(new_point.x,new_point.y,'gx')
                #self.ID_jet_boundary.append(self.ID)
                self.ID +=1
                self.contour_fan = 0   
                self.add_break_ID()
                self.new_fan = 1
                # if (self.centre_line_intercept):
                #     self.END_SIM = 1 
            else:
                if(self.chr_array[ID_temp].pt_type!="same_fam"):
                    temp1 = self.same_fam_point(self.chr_array[self.ID_right_chr[0]],self.chr_array[ID_temp])
                    temp2 = self.general_point(self.chr_array[ID_temp],self.chr_array[self.ID-1]) 
                    if (temp1.x < temp2.x) and (temp1.x>self.chr_array[ID_temp].x) :
                        #self.ID_right_chr.pop(-1); self.ID_left_chr.pop(-1)
                        #self.compression_offset += 1
                        #self.plot_chr=1
                        self.ID_left_chr.pop(-1)
                        self.ID_right_chr.pop(-1)
                        #print(temp1.x)
                        if (self.plot_chr):
                            plt.plot(self.chr_array[self.ID_right_chr[0]].x,self.chr_array[self.ID_right_chr[0]].y,'bx',self.chr_array[ID_temp].x,self.chr_array[ID_temp].y,'rx')
                            plt.plot(temp1.x,temp1.y,'go')

                            #plt.plot([self.chr_array[ID_temp].x, temp1.x],[self.chr_array[ID_temp].y, temp1.y])
                        self.chr_array = np.append(self.chr_array,temp1)
                  
                        new_point = self.general_point(self.chr_array[-1],self.chr_array[self.ID-1],plot_chr=self.plot_chr) 
                        #plt.plot(self.chr_array[self.ID-1].x,self.chr_array[self.ID-1].y,'ro')
                        self.ID += 1      
                    else:
                        new_point = temp2
                        if (new_point.x<=self.spike.length):
                            pass
                        if (self.plot_chr):
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

        ## END OF MOC SECTION #####
        # function order is important
        if(clean_mesh):
            self.clean_data()
        self.to_arrays()
        self.calc_flow_properties()
        #self.compute_thrust()

    def compute_initial_expansion_fan(self):
        M_max = gd.PR_expansion_mach(self.PR,self.gamma)
        # print('M_max: ' + str(M_max))
        mach_fan = np.linspace(1.1,M_max,self.n)

        (T_ratio,p_ratio,rho_ratio,a_ratio) = gd.isentropic_ratios(0,mach_fan,self.gamma)

        V_fan = a_ratio*self.spike.a_c*mach_fan

        W_fan = V_fan/self.V_l

        theta_fan = -gd.prandtl_meyer(mach_fan,self.gamma) + self.slope_init

        angle_fan = gd.mach_angle(mach_fan)

        # print(180/np.pi*np.arcsin(np.sqrt((gamma-1)/2*(1/W_fan**2-1))))
        # print(W_fan)

        # print(angle_fan*180/np.pi)

        x_fan = np.ones(angle_fan.shape)*self.spike.lip_x

        y_fan = np.ones(angle_fan.shape)*self.spike.lip_y

        #print(theta_fan*180/np.pi)
        # print(gd.mach_angle_velocity_ratio(gd.prandtl_meyer(2.3,gamma),0.3,gamma))       

        initial_point = self.contour_point(chr_point(self.gamma,x_fan[0],y_fan[0],theta_fan[0],W_fan[0],'N/A'),plot_chr=self.plot_chr)
        self.ID += 1
        self.ID_contour_chr.pop(0)
        self.chr_array = np.append(self.chr_array,initial_point)

        for point in x_fan[1:-1]:
            temp_point = chr_point(self.gamma,x_fan[self.ID],y_fan[self.ID],theta_fan[self.ID],W_fan[self.ID],'N/A')
            new_point = self.general_point(temp_point,self.chr_array[self.ID-1],plot_chr=self.plot_chr)
            # adding to arrays
            self.chr_array = np.append(self.chr_array,new_point)
            self.ID += 1

        first_jet = chr_point(self.gamma,x_fan[-1],y_fan[-1],theta_fan[-1],W_fan[-1],'N/A')
        second_jet = self.jet_boundary_point(first_jet,self.chr_array[self.ID-1],plot_chr=self.plot_chr)
        self.chr_array = np.append(self.chr_array,second_jet)
        #self.ID_jet_boundary.append(self.ID)
        self.ID += 1
        self.add_break_ID()

    def general_point(self,A,B,plot_chr=0):
        # given points A & B (B.y>A.y) in char mesh, computes 3rd point

        x_c = (A.x*np.tan(A.theta+A.mu)-A.y+B.y-B.x*np.tan(B.theta-B.mu))/(np.tan(A.theta+A.mu)-np.tan(B.theta-B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta+A.mu) + A.y
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta-B.mu)

        if not (self.on_nozzle_contour(A)):
            theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A*(x_c-A.x)/A.y)+B.W+2*B.W*B.theta*np.tan(B.mu))/(A.W*np.tan(A.mu)+2*B.W*np.tan(B.mu))
        else:
            theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A/A.y*(x_c-A.x))+B.W+B.W*(B.theta*np.tan(B.mu)+m_B/B.y*(x_c-B.x)))/(A.W*np.tan(A.mu)+B.W*np.tan(B.mu))
        W_c = A.W + A.W*(np.tan(A.mu*(theta_c-A.theta))+l_A/A.y*(x_c-A.x))

        # checking if point is greater than length
        self.ID_right_chr.append(self.ID)  
        self.ID_left_chr.append(self.ID)        
            
        
        if (x_c <= self.spike.length):
            pass
            # self.ID_left_chr.append(self.ID)
            # self.ID_right_chr.append(self.ID)              
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'k',[B.x,x_c],[B.y,y_c],'k')
      
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c,'general')

    def same_fam_point(self,A,B,plot_chr=0):
        # DEPRECATED, SHOULD BE CHECKED
        # given points A & B (B.y>A.y) in char mesh, computes 3rd point in the case that the interception is of the same family of points
        with np.errstate(invalid='ignore',divide='ignore'):
            x_c = (A.x*np.tan(A.theta+A.mu)-A.y+B.y-B.x*np.tan(B.theta+B.mu))/(np.tan(A.theta+A.mu)-np.tan(B.theta+B.mu))
            y_c = (x_c-A.x)*np.tan(A.theta+A.mu) + A.y
            l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
            m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta+B.mu)
            theta_c = (-A.W-A.W*(-A.theta*np.tan(A.mu)+l_A/A.y*(x_c-A.x))+B.W+B.W*(B.theta*np.tan(B.mu)+m_B/B.y*(x_c-B.x)))/(A.W*np.tan(A.mu)+B.W*np.tan(B.mu))
            W_c = A.W + A.W*(np.tan(A.mu*(theta_c-A.theta))+l_A/A.y*(x_c-A.x))
        # if (self.plot_chr):
        #     plt.plot([A.x,x_c],[A.y,y_c],'r',[B.x,x_c],[B.y,y_c],'r')
        if (x_c <= self.spike.length):
            pass
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'k',[B.x,x_c],[B.y,y_c],'k')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c,'same_fam')

    def contour_point(self,A,plot_chr=0):
        # Given a chr_point A, computes a chr_point that intersects with the nozzle boundary
        line = lambda x: np.tan(A.theta+A.mu)*(x-A.x) + A.y
        # print('line: ' + str((A.theta+A.mu)*180/np.pi))
        intercept = lambda x: interpolate.splev(x,self.tck,der=0) - line(x)
        x_c = optimize.brentq(intercept,self.spike.x[0],self.spike.x[-1])
        y_c = line(x_c)  
        theta_c = np.arctan(interpolate.splev(x_c,self.tck,der=1)) 
        l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)  
        W_c = A.W + A.W*((theta_c-A.theta)*np.tan(A.mu)+l_A/A.y*(x_c-A.x))
        #plt.plot(x_c,y_c,'rX')
        # self.ID_contour_chr.append(self.ID)
        self.ID_contour_chr.append(self.ID)  
        if (x_c <= self.spike.length):  
            pass  
        if(plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'k')

        return chr_point(self.gamma,x_c,y_c,theta_c,W_c,'contour')

    def jet_boundary_point(self,A,B,plot_chr=0):
        # given points A and B (where A is the previous jet boundary point and A.x < B.x) return the next boundary point
        x_c = (A.x*np.tan(A.theta)-A.y+B.y-B.x*np.tan(B.theta-B.mu))/(np.tan(A.theta)-np.tan(B.theta-B.mu))
        y_c = (x_c-A.x)*np.tan(A.theta) + A.y
        mu_c = A.mu
        W_c = A.W
        m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta-B.mu)
        theta_c = B.theta + (-(W_c-B.W)/B.W + m_B*(x_c-B.x)/B.y)/np.tan(B.mu)
        
   
        self.ID_jet_boundary.append(self.ID)

        if (x_c <= self.spike.length):
            pass
            # self.ID_jet_boundary.append(self.ID)
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'r',[B.x,x_c],[B.y,y_c],'r')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c,'jet')
    def internal_jet_point(self,A,B,plot_chr=0):
        # given points A and B (where A is previous internal jet boundary point)
        x_c = (A.x*np.tan(A.theta)-A.y+B.y - B.x*np.tan(B.theta + B.mu))/(np.tan(A.theta) - np.tan(B.theta + B.mu))
        y_c = (x_c - A.x)*np.tan(A.theta) + A.y
        if (y_c < 0):
            W_c = A.W
            mu_c = A.mu 
            m_B = np.sin(B.mu)*np.sin(B.theta)*np.tan(B.mu)/np.cos(B.theta-B.mu)
            theta_c = B.theta + (-(W_c -B.W)/B.W + m_B*(x_c - B.x)/B.y)/np.tan(B.mu)
        else:
            #centre point
            self.centre_line_intercept = 1
            x_c = A.x - A.y/(np.tan(A.theta))   # + A.theta
            y_c = 0  
            theta_c = 0   
            W_c = A.W
            # print('W_c = ' + str(W_c))
            # x_c = A.x - A.y/(np.tan(A.theta+ A.theta))
            # y_c = 0  
            # theta_c = 0   
            # #print('W_c = ' + str(W_c))
            # l_A = np.sin(A.mu)*np.sin(A.theta)*np.tan(A.mu)/np.cos(A.theta+A.mu)
            # W_c = A.W + A.W*(np.tan(A.mu*(theta_c-A.theta))+l_A/A.y*(x_c-A.x))

        
        self.ID_contour_chr.append(self.ID)
        if (x_c <= 2*self.spike.length):
            pass
        if (plot_chr):
            plt.plot([A.x,x_c],[A.y,y_c],'g',[B.x,x_c],[B.y,y_c],'g')
        return chr_point(self.gamma,x_c,y_c,theta_c,W_c,'contour')                

    def add_break_ID(self):
        #self.ID_left_chr.append(-1)
        self.ID_right_chr.append(-1)
        #self.ID_jet_boundary.append(-1)
        #self.ID_contour_chr.append(-1)      

    def base_pressure(self):
        ###
        # UNIVERSITY OF ROME MODEL (2002)
        ###
        phi = np.arctan(interpolate.splev(self.spike.x[0], self.tck, der=1))
        PHI = (-0.2*phi**4-5.89*phi**2+20179.84)/(phi**4+20179.84)

        return self.p_atm*(0.05 + 0.967)**PHI

    def on_nozzle_contour(self,chr_point_obj):
        line = lambda x: np.tan(chr_point_obj.theta+chr_point_obj.mu)*(x-chr_point_obj.x) + chr_point_obj.y
        intercept = lambda x: interpolate.splev(x,self.tck,der=0) - line(x)
        return np.sign(intercept(self.spike.x[0])*intercept(self.spike.x[-1])) <= 0

    def chr_point_less_zero(self):
        if len(self.ID_contour_chr) > 0:
            return self.chr_array[self.ID_contour_chr[-1]].y < 0
        else:
            return 1

    def within_x_bound(self):
        if len(self.ID_contour_chr) > 0:
            return self.chr_array[self.ID_contour_chr[-1]].x < spike.x.max()*20
        else:
            return 1

    def contour_converge(self):
        if len(self.ID_contour_chr) > 1:
            return np.abs(self.chr_array[self.ID_contour_chr[-1]].y) < np.abs(self.chr_array[self.ID_contour_chr[-2]].y)
        else:
            return 1
    def flip_plug(self):
        # confirms that the plug y values are negative and flips if required
        if self.spike.lip_y > 0:
            self.spike.lip_y = self.spike.lip_y*-1
            self.spike.y = self.spike.y*-1



    def to_arrays(self):
        self.x = np.array([]); self.y = np.array([]); self.M = np.array([]); self.mu = np.array([]); self.V = np.array([]); self.theta = np.array([])
        for point in self.chr_array:
            #print(type(point.x))
            self.x = np.append(self.x,point.x)
            self.y = np.append(self.y,point.y)
            self.M = np.append(self.M,point.M)
            self.mu = np.append(self.mu,point.mu)
            self.V = np.append(self.V,point.W*self.V_l)
            self.theta = np.append(self.theta,point.theta)

        #takes chr_array obj points to arrays
    def calc_flow_properties(self):
        T_ratio,p_ratio,rho_ratio,a_ratio = gd.isentropic_ratios(0,self.M,self.gamma)
        if(p_ratio > 1):
            print('ERROR')
        self.T = self.spike.T_c*T_ratio
        self.p = self.spike.p_c*p_ratio   
        self.a = self.spike.a_c*a_ratio
        self.rho = self.spike.rho_c*rho_ratio


    def clean_data(self):

        def clean_ID_list(chr_array):
            contour_ID= []
            jet_ID = []
            ID = 0
            for point in chr_array:
                if point.pt_type == 'contour':
                    contour_ID.append(ID)
                elif point.pt_type == 'jet':
                    jet_ID.append(ID)

                ID += 1

            return contour_ID, jet_ID 

        curr_ID = 0
        del_ID = []

        # creating list of IDs of all points to be deleted from mesh
        for point in self.chr_array:
            if point.x < self.chr_array[0].x or point.y > 0 or point.x > self.chr_array[self.ID_contour_chr[-1]].x*self.downstream_factor or np.isnan(point.x):
                del_ID.append(curr_ID)
            curr_ID += 1


        # deleting IDs from lists
        self.chr_array = np.delete(self.chr_array,del_ID)
        

        #print(len(self.ID_jet_boundary))
        self.ID_contour_chr,self.ID_jet_boundary = clean_ID_list(self.chr_array)

    def compute_thrust(self,approx_method,n):       
        # # constructing spline representation for jet boundary
        jet_bound_x = np.concatenate((np.array([self.spike.lip_x]),self.x[self.ID_jet_boundary]))
        jet_bound_y = np.concatenate((np.array([self.spike.lip_y]),self.y[self.ID_jet_boundary]))
        
        # filter for uniques
        jet_bound_x, indices = np.unique(jet_bound_x, return_index=True)
        jet_bound_y = jet_bound_y[indices]

        try:
        #constructing jet boundary spline
            tck_jet_bound = interpolate.splrep(jet_bound_x,jet_bound_y)
        except TypeError:
            print(self.altitude)
            print(str(jet_bound_x))
            print(str(jet_bound_y))
       
        # constructing plane on which to evaluate expanded gas properties
        x_plane = np.ones(n,)*self.x[self.ID_contour_chr[-1]];
        y_points = np.linspace(0,interpolate.splev(self.x[self.ID_contour_chr[-1]],tck_jet_bound),n);
        #plt.plot(x_plane,y_points,'ro')
        # constructing rbf functions for interpolation of properties
        V_grid = interpolate.griddata((self.x,self.y),self.V,(x_plane,y_points),method=approx_method) # nearest and cubic may also work 
        T_grid = interpolate.griddata((self.x,self.y),self.T,(x_plane,y_points),method=approx_method)
        theta_grid = interpolate.griddata((self.x,self.y),self.theta,(x_plane,y_points),method=approx_method)
        rho_grid = interpolate.griddata((self.x,self.y),self.rho,(x_plane,y_points),method=approx_method)
        P_grid = interpolate.griddata((self.x,self.y),self.p,(x_plane,y_points),method=approx_method)
        #computing thrust
        Ve_grid = V_grid*np.cos(theta_grid) # V*cos(theta)*r *dr 
        # fig1, (ax1) = plt.subplots(1,1)
        # vel_plot = ax1.scatter(y_points,Ve_grid)
        # vel_label = 'Ve ' + str(self.altitude)
        # print(vel_label)
        # plt.plot(y_points,Ve_grid,label = vel_label)

        tck_contour = interpolate.splrep(self.spike.x,self.spike.y)

        

        on_nozzle_idx = self.x[self.ID_contour_chr] < self.spike.x.max()
        on_contour_x = self.x[self.ID_contour_chr][on_nozzle_idx]
        on_contour_y = self.y[self.ID_contour_chr][on_nozzle_idx]


        # print('thurst!')
        # plt.plot(on_contour_x,on_contour_y,'x')
        # plt.show()
        contour_der = interpolate.splev(on_contour_x,tck_contour,der=1)

        contour_angle = np.arctan(contour_der)
        

        P_contour = interpolate.griddata((self.x,self.y),self.p,(on_contour_x,on_contour_y),method=approx_method) - self.p_atm

        P_2D = P_contour*np.sqrt(contour_der**2+1)*2*np.pi*np.sin(contour_angle)*2*np.pi*on_contour_y
        # p_label = 'P ' + str(self.altitude) + ', patm:' + str(self.p_atm)
        # plt.plot(y_points,np.cos(theta_grid),label =p_label)
        # plt.legend()
        #plt.show()
        A = y_points*2*np.pi
        thrust_grid = rho_grid*Ve_grid**2*A ## CHECK THIS!!!!!!!
        thrust_momentum = np.trapz(thrust_grid,y_points) # check with emerson
        #thrust_pressure = np.trapz((P_2D),on_contour_x) # check with emerson

        return thrust_momentum #+ thrust_pressure
        #FUN PLOTS
        # print(thrust_momentum)
        # print(thrust_pressure)
        # fig1, (ax1,ax2) = plt.subplots(1,2)
        # vel_plot = ax1.scatter(self.x,self.y,c=(self.V*np.cos(self.theta)),cmap=cm.coolwarm)
        # ax1.plot(x_plane[-1],y_points[-1],'x')
        # vel_interp = ax2.plot(y_points,Ve_grid)
        # ax1.axis('equal')
        # ax2.set_xlim(y_points.min(),y_points.max())
        # plt.colorbar(vel_plot,ax = ax1)


# #### NOZZLE INITIAL CONDITIONS
# r_e = 0.067/2 #0.034 # likely too large
# expansion_ratio = 6.64 #8.1273
# A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
# gamma = 1.2381 #np.mean([1.2534,1.2852])
# T_c = 2833.63
# p_c = 34.474*10**5
# rho_c = 3.3826
# R = (1-1/gamma)*1.8292*1000#1.8292 1.6196
# a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 
# n = 1000

# ### DESIGN OF SPIKE

# spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 1.0)

# spike.y = spike.y*-1
# spike.lip_y = spike.lip_y*-1

# line = lambda x: np.tan(-np.pi/3)*(x-spike.lip_x) + spike.lip_y


# ## METHOD FOR CALCULATING FUNCTION INTERCEPTS

# tck = interpolate.splrep(spike.x,spike.y,full_output=0)


# MOC_mesh = chr_mesh(spike,gamma,0,50,downstream_factor=1.2,plot_chr=1)

# #plt.plot(MOC_mesh.x,MOC_mesh.y,'ro')

# #plt.plot(MOC_mesh.x[MOC_mesh.ID_jet_boundary],MOC_mesh.y[MOC_mesh.ID_jet_boundary],'go-')
# MOC_mesh.compute_thrust('nearest',10)
# #plt.plot(MOC_mesh.x[MOC_mesh.ID_contour_chr],MOC_mesh.y[MOC_mesh.ID_contour_chr],'bo-')

# #plt.plot(MOC_mesh.x[MOC_mesh.ID_contour_chr[-1]],0,'bo')
# #plt.plot(spike.x,interpolate.splev(spike.x,tck,der=0),spike.lip_x,spike.lip_y,'X')
# #plt.axis('equal')
# # fig = plt.figure()
# # ax = fig.add_subplot(111)
# # ax.scatter(MOC_mesh.x,MOC_mesh.y,c=MOC_mesh.M, cmap = cm.coolwarm)
# # #plt.plot(MOC_mesh.x,MOC_mesh.y,'r.')
# # # cax = ax.imshow(MOC_mesh,interpolation='nearest',cmap=cm.coolwarm)

# # # cbar = fig.colorbar(cax)
# # plt.colorbar(ax=ax)
# #plt.axis('equal')
# plt.show()


# M = optimize.brentq(lambda M: gd.expansion_ratio_zero(1,M,gamma,expansion_ratio),1,10)
