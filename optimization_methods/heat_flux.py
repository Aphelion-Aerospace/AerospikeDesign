#Uses Mayer's method to calculate heat convection coefficient, then calculate total heat flux
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
#from design_nozzle import spike
from math import pi

#input variables from NASA CEA in metric units:
Pr=0.55645 #average throat to exit Prandtl's number
Cp=1.724 #[KJ/KG-K] average throat to exit constant pressure heat capacity
Gamma=1.2689 #average throat to exit ratio of specific heat
c=0.003883468 #[millipoise/K^w] viscocity to temperature coefficient
w=0.678083301 #viscocity to temperature exponent
To=2833 #[K] stagnation/chamber temperature

#user input variable in metric units:
T_w=700 #[K] desired temperature of nozzle wall

def heat_flux(Pr,Cp,Gamma,c,w,To,spike):
    class heat_transfer:
        def __init__(self,pr,cp,gamma,c,w,To,T_w,s,r,M,u,rho,T_inf):
            #input parameters and unit conversion if needed
            self.pr=pr
            self.cp=cp*0.238845 #[Btu/lb-R]
            self.gamma=gamma
            self.c=c*0.0000671968994813/1.8**w #[lb/ft-s-R^w]
            self.w=w
            self.To=To*1.8 #[R]
            self.T_w=T_w*1.8 #[R]
            self.s=s*3.28084 #[ft] nozzle wall length
            self.r=r*3.28084 #[ft] nozzle radius
            self.M=M #Mach number
            self.u=u*3.28084 #[ft] free stream velocity
            self.rho=rho*0.062428 #[lb/ft^3] free stream density
            self.T_inf=T_inf*1.8 #[R] free stream temperature

            #calculate additonal parameters
            self.mu=self.c*(self.T_inf**self.w) #[lb/ft-s] free stream viscocity
            self.T_ad=self.To*(1+(self.gamma-1)/2*self.M**2*self.pr**(1/3))/(1+(self.gamma-1)/2*self.M**2) #[R] adiabatic temperature
            self.T_star=0.5*(self.T_w+self.T_inf)+0.22*(self.T_ad-self.T_inf) #[R]
            self.b=(self.T_inf/self.T_star)**(0.8-0.2*self.w) # beta factor

            #initiate arrays to calculate h2 and h3
            self.h2_den=self.s*0
            self.h3_den=self.s*0

            #intergrate to get the denominator of h2 and compute h2
            self.h2_y=self.b**(5/4)*self.rho*self.u/self.mu
            i=1
            while i <= len(self.s):
                self.h2_den[i-1]=np.trapz(self.h2_y[:i],self.s[:i])
                i=i+1
            self.h2_den[0]=1 #to prevent division by 0
            self.h2=0.0296*self.b**(5/4)*self.pr**(-2/3)*self.rho*self.cp*self.u/self.h2_den**0.2

            #intergrate to get denominator of h3 and calculate h3
            self.h3_y=(self.r*self.b)**(5/4)*self.rho*self.u/self.mu
            i=1
            while i <= len(self.s):
                self.h3_den[i-1]=np.trapz(self.h3_y[:i],self.s[:i])
                i=i+1
            self.h3_den[0]=1 #to prevent division by 0
            self.h3=self.h2*(self.r**(5/4)*self.h2_den/self.h3_den)**0.2 #[BTU/ft^2-sec-R] convective heat transfer coefficient
            self.h3[0]=0 #the first h3 will be infinity, so input 0 in instead
            self.h3_in=self.h3/144 #[BTU/in^2-sec-R] unit conversion to in^2
            self.h3_m=self.h3*20441.748028012 #[W/m^2-K] unit conversion to metric

    heat=heat_transfer(Pr,Cp,Gamma,c,w,To,T_w,spike.s,spike.y,spike.M,spike.V,spike.rho,spike.T)

    #compute heat flux assuming cylinder with height dx
    i=277 #skip the first 0.4mm because h3 approached infinity
    q=np.zeros(len(spike.x)-i)
    q[0]=heat.h3_m[i]*(heat.T_ad[i]/1.8-T_w)*spike.x[i]*(2*pi*spike.y[i])
    i=i+1
    while i < len(spike.x):
        q[i-277]=heat.h3_m[i]*(heat.T_ad[i]/1.8-T_w)*(spike.x[i]-spike.x[i-1])*(2*pi*spike.y[i])
        i=i+1

    #compute total heat flux
    return total_q=np.sum(q) #[W]

# print ("Total heat flux ="+ str(heat_flux(Pr,Cp,Gamma,c,w,To))+ "Watt")

# # csv_array = np.array([heat.h3_in,heat.h3_m]) #export heat transfer coefficient to excel
# # np.savetxt('heat_transfer_coefficient.csv', csv_array.T, delimiter = ',')

# plt.plot(spike.x[1:]*100,heat.h3_m[1:])
# plt.ylabel('Heat Transfer Coefficient (W/m^2-K)')
# plt.xlabel('Nozzle length (cm)')
# plt.show()
