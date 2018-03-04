import numpy as np 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import scipy.optimize
from scipy import interpolate
import multiprocessing as mp
import copy
import pickle
from matplotlib import cm 

import gasdynamics as gd
from heat_flux import heat_flux
from plug_nozzle_angelino import plug_nozzle 
import MOC
import MOC_its


## NASA CEA CONSTANTS
class CEA_constants():
	def __init__(self,altitude):
		"""Class that should eventually queery NASA CEA data, for now holds all values constant

		Args:
			altitude (float): altitude in meters

		Returns:
			Stores NASA CEA data in atributes
		"""	
		self.gamma = 1.2381 #np.mean([1.2534,1.2852])
		self.T_c = 2833.63 # combustion chamber temperature
		self.p_c = 34.474*10**5 # combustion chamber pressure
		self.rho_c = 3.3826 # combustion chamber density 
		self.a_c = np.sqrt(self.gamma*(1-1/self.gamma)*200.07*self.T_c) # combustion chamber sound speed
		self.Pr = 0.55645 #average throat to exit Prandtl's number
		self.cp = 1.724 #[KJ/KG-K] average throat to exit constant pressure heat capacity
		self.c = 0.003883468 #[millipoise/K^w] viscocity to temperature coefficient
		self.w = 0.678083301 #viscocity to temperature exponent


class aerospike_optimizer():
	"""Class intended to be used to optimize an aerospike nozzle

	Args:
		r_e (float): outer radius of aerospike nozzle
		T_w (float): desired nozzle temperatures
		alpha (float): weighting of thrust
		beta (float): weighting of cooling requirement
		design_alt_init (float): initial design altitude
		truncate_ratio_init (float): inital truncation length
		chr_mesh_n (int): number of expansion waves for MOC
		no_alt_range (float): number of altitude ranges to optimize over
		no_core (int>0): number of cores to be used in computation

	Returns:
		cost_opt_contour_params (float): cost computed with spike.x and spike.y as input parameters
		cost_opt_design_params (float): cost computed with design_alt and truncate_ratio as input parameters
	"""
	def __init__(self,r_e,T_w,alpha,beta,design_alt_init,truncate_ratio_init,chr_mesh_n=120,no_alt_range = 30,no_core=1):
		self.r_e = r_e
		self.T_w = T_w
		self.alpha = alpha
		self.beta = beta
		self.design_alt_init = design_alt_init
		self.truncate_ratio_init = truncate_ratio_init
		self.chr_mesh_n = chr_mesh_n
		self.no_alt_range_int = int(no_core*round(no_alt_range/no_core))
		self.no_core = no_core

		self.spike_init = self.__design_angelino_nozzle(design_alt,truncate_ratio,CEA,r_e)

		self.spike_opt = copy.deepcopy(self.spike_init)

		self.CEA = CEA_constants(0)

	def __compute_thrust_over_range(self,plug_nozzle_class,alt_range,gamma,send_end,downstream_factor=1.2,chr_mesh_n=50):
		"""Function that computes the thrust at each point over a specified range, intended to be used as apart of a process

		Args: 
			plug_nozzle_class (class plug_nozzle_angelino.plug_nozzle): the plug nozzle for which the thrust is to be computed
			alt_range (np.array): discrete altitudes over which thrust should be computed
			gamma (float): ratio of specifc heats
			send_end (class multiprocessing.connection.Connection): part of recv_end, send_end pair indicating which recv_end to send to
			downstream_factor (float): specifies what factor downstream of the final contour point MOC should be continued up to
			chr_mesh_n (int>0): specifies number expansion waves at lip for MOC

		Returns:
			None
			send_end.send(thrust_range); thrust_range (np.array): contains thrust computed at each alt_range point

		"""
		thrust_range = np.zeros(alt_range.shape)
		for i in range(alt_range.shape[0]):
			# print(alt_range[i])
			try:
				MOC_mesh = MOC.chr_mesh(plug_nozzle_class,gamma,alt_range[i],chr_mesh_n,downstream_factor=downstream_factor)
				thrust_range[i] = MOC_mesh.compute_thrust('nearest',30)
				print(thrust_range[i])
			except:
				thrust_range[i] = 0
		send_end.send(thrust_range)

	def __multicore_thrust_compute(self,plug_nozzle_class,altitude_range,gamma,downstream_factor=1.2,chr_mesh_n=50,no_core=1):
		proc_list = []
		pipe_list =[]

		alt_range_split = np.split(altitude_range,no_core)

		for i in range(no_core):
			recv_end, send_end = mp.Pipe(False)
			args = (plug_nozzle_class,alt_range_split[i],gamma,send_end,downstream_factor,chr_mesh_n)
			proc = mp.Process(target=self.__compute_thrust_over_range, args = args)
			proc_list.append(proc)
			pipe_list.append(recv_end)
			proc.start()

		for proc in proc_list:
			proc.join()

		thrust_range = [x.recv() for x in pipe_list]

		thrust_range = np.concatenate(thrust_range)

		# for thread in threads:
		# 	thread.map

		return thrust_range

	def __design_angelino_nozzle(self,design_alt,truncate_ratio,CEA,r_e):
		(p_atm,T_atm,rho_atm) = gd.standard_atmosphere([design_alt])

		PR = CEA.p_c/p_atm

		M_e = gd.PR_expansion_mach(PR,CEA.gamma)

		expansion_ratio = gd.expansion_ratio(1,M_e,CEA.gamma)#6.64 #8.1273
		# print('Exp. ratio: ' + str(expansion_ratio))
		# print('PR: ' + str(PR))

		A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)

		return plug_nozzle(expansion_ratio,A_t,r_e,CEA.gamma,CEA.T_c,CEA.p_c,CEA.a_c,CEA.rho_c,100,truncate_ratio = truncate_ratio)

	def __cost_end_func(self,no_alt_range_int,spike,CEA,downstream_factor=1.2,chr_mesh_n=120,no_core=1):
		alt_range = np.linspace(0,9144,no_alt_range_int)

		# shuffle arrays so each core computes similar complexity on average
		np.random.shuffle(alt_range)
		(p_atm_r,T_atm_r,rho_atm_r) = gd.standard_atmosphere(alt_range)

		thrust_range = self.__multicore_thrust_compute(spike,alt_range,CEA.gamma,downstream_factor=1.2,chr_mesh_n=chr_mesh_n,no_core=no_core)

		# unshuffle arrays
		ordered_idx = np.argsort(alt_range)
		alt_range = alt_range[ordered_idx]; thrust_range = thrust_range[ordered_idx]

		work = np.trapz(thrust_range,alt_range)
		# plt.plot(alt_range,thrust_range,'o')
		# plt.show()
		print('work = ' + str(work))

		## heat transfer required
		total_heat_flux = heat_flux(CEA.Pr,CEA.cp,CEA.gamma,CEA.c,CEA.w,CEA.T_c,T_w,spike)

		return (work, total_heat_flux)

	def __cost_func_contour_params(self,params,spike,T_w,CEA,alpha,beta,chr_mesh_n,no_alt_range_int,no_core):
		params = np.asarray(params)
		if len(params.shape) > 1:
			raise ValueError('Input params not correct shape. Should be 1D array')
		x_vals,y_vals = np.split(params,2)	


		spike.x = x_vals; spike.y = y_vals


		work, total_heat_flux = self.__cost_end_func(no_alt_range_int,spike,CEA,downstream_factor=1.2,chr_mesh_n=chr_mesh_n,no_core=no_core)
		print(total_heat_flux)
		return -alpha*work + beta*total_heat_flux
		
	def __cost_func_design_params(self,params,T_w,CEA,r_e,alpha,beta,chr_mesh_n,no_alt_range_int,no_core):
		params = np.asarray(params)
		if len(params.shape) > 1:
			raise ValueError('Input params not correct shape. Should be 1D array')
		design_alt, truncate_ratio = np.split(params,2)
		spike = self.__design_angelino_nozzle(design_alt,truncate_ratio,CEA,r_e)

		work, total_heat_flux = self.__cost_end_func(no_alt_range_int,spike,CEA,downstream_factor=1.2,chr_mesh_n=chr_mesh_n,no_core=no_core)
		print('work = ' + str(work))
		return -alpha*work + beta*total_heat_flux


	def multicore_animate_to_its(self,max_its,no_core):
		iterations = np.arange(max_its)
		np.random.shuffle(iterations)
		iterations_split = np.split(iterations,no_core)
		for i in range(no_core):
			args = (iterations_split[i],)
			proc = mp.Process(target=self.animate_to_its, args = args)
			proc.start()



	def animate_to_its(self,num_its):
		spike = self.__design_angelino_nozzle(self.design_alt_init,self.truncate_ratio_init,self.CEA,self.r_e)


		for its in num_its:
			spike.plot_contour(plt)
			spike.y = spike.y*-1; spike.lip_y = spike.lip_y*-1
			spike.plot_contour(plt)

			try:
				MOC_mesh = MOC_its.chr_mesh(spike,self.CEA.gamma,self.design_alt_init,30,its,downstream_factor=1.1,plot_chr=1)
		
			except:
				pass

			plt.ylim((-0.1,0.1))
			plt.xlim((-0.01,0.26))
			name = 'animation_1_full_length/fig' + str(its)
			plt.savefig(name)
			plt.close()

	def animate_over_range(self,plug_nozzle_class,CEA,alt_range,flight_range,downstream_factor=1.2,chr_mesh_n=50):
		
		for i in range(alt_range.shape[0]):

			fig, (ax1,ax2) = plt.subplots(1,2)#,gridspec_kw={'width_ratios':[8,1]})
			MOC_mesh = MOC.chr_mesh(plug_nozzle_class,CEA.gamma,alt_range[i],chr_mesh_n,downstream_factor=5,plot_chr=0,clean_mesh=1)
		
			plug_nozzle_class.plot_contour(ax1)
			plug_nozzle_class.y = plug_nozzle_class.y*-1; plug_nozzle_class.lip_y = plug_nozzle_class.lip_y*-1
			plug_nozzle_class.plot_contour(ax1)
			ax1.plot([plug_nozzle_class.x[-1],plug_nozzle_class.x[-1]],[plug_nozzle_class.y[-1],plug_nozzle_class.y[-1]*-1],'b-')
			ax1.set_ylim([-0.1,0.1])#ax1.set_ylim([-0.05,0.05])# 
			ax1.set_xlim([-0.01,0.26])#ax1.set_xlim([-0.01,0.05])# 

			contourf_grid = 1000

			x_plt = np.linspace(MOC_mesh.x.min(),MOC_mesh.x.max(),contourf_grid)
			y_plt = np.linspace(MOC_mesh.y.min(),MOC_mesh.y.max(),contourf_grid)
			X_plt,Y_plt = np.meshgrid(x_plt,y_plt)
		
			M_contour=interpolate.griddata((MOC_mesh.x,MOC_mesh.y),MOC_mesh.T,(X_plt,Y_plt),method='linear')
			inner_contour_tck = interpolate.splrep(MOC_mesh.x[MOC_mesh.ID_contour_chr],MOC_mesh.y[MOC_mesh.ID_contour_chr])
			def contour_curve(x):
				return interpolate.splev(x,inner_contour_tck)
			M_contour[Y_plt > contour_curve(X_plt)] = np.nan 

			M_fill = ax1.contourf(X_plt,Y_plt,M_contour,cmap=cm.jet)

			#v=np.linspace(0,5,10)
			#ax1.contourf(X_plt,-1*Y_plt,M_contour,cmap=cm.jet)
			plt.colorbar(M_fill,ax=ax1)
			#plt.clim(0,5)
			name = 'animation_4_full_length_p/fig' + str(int(alt_range[i]))
			#title = str(int(alt_range[i]*3.28084)) + str('ft')
			#plt.suptitle(title,fontsize=40)
			ax1.set_aspect('equal','box')

			#ax2.xaxis.set_visible(False); 

			x_flight = np.linspace(0,flight_range[i])

			# ax2.plot(plug_nozzle_class.x*1000,plug_nozzle_class.T)
			# ax2.set_xlabel('x (mm)')
			# ax2.set_ylabel('Temperature (K)')
			# ax2.plot(x_flight,60000*self.__fake_flight_parabola(x_flight))
			# ax2.plot(flight_range[i],60000*self.__fake_flight_parabola(flight_range[i]),'rx')
			# ax2.set_ylim([0,60000])
			# ax2.set_yticks([0,10000,20000,30000,40000,50000,60000])
			# ax2.set_xlim([0,1])
			fig.set_size_inches(18.5,10.5)
			plt.savefig(name,dpi=100)

			plt.close()	


	def __fake_flight_parabola(self,x):
		return -(x-1)**2 + 1	

	def multicore_animate_range(self,downstream_factor=1.2,chr_mesh_n=50,no_core=1):
 #3.28084

		alt_range = np.linspace(0,16764,self.no_alt_range_int)
		flight_range = np.linspace(0,1,self.no_alt_range_int)
		rand_idx = np.random.permutation(len(alt_range))
		alt_range = alt_range[rand_idx]; flight_range = flight_range[rand_idx]

		alt_range_split = np.split(alt_range,no_core)
		flight_range_split = np.split(flight_range,no_core)
		plug_nozzle_class = self.__design_angelino_nozzle(self.design_alt_init,self.truncate_ratio_init,self.CEA,self.r_e)
		

		for i in range(no_core):

			args = (plug_nozzle_class,self.CEA,alt_range_split[i],flight_range_split[i],downstream_factor,chr_mesh_n)
			proc = mp.Process(target=self.animate_over_range, args = args)
			proc.start()



	def cost_opt_contour_params(self,params): 
		"""params = np.concatenate((spike_opt.x,spike_opt.y))
		"""
		return self.__cost_func_contour_params(params,self.spike_opt,self.T_w,self.CEA,self.alpha,self.beta,self.chr_mesh_n,self.no_alt_range_int,self.no_core)

	def cost_opt_design_params(self,params): 
		"""params = [design_alt,truncation_ratio]
		"""
		return self.__cost_func_design_params(params,self.T_w,self.CEA,self.r_e,self.alpha,self.beta,self.chr_mesh_n,self.no_alt_range_int,self.no_core)
if __name__ == '__main__':

	## CONSTANTS OF DESIGN FOR AERODYNAMICS
	r_e = 0.027 #0.034 # likely too large

	## CONSTANTS OF DESIGN FOR HEAT FLUX
	#user input variable in metric units:
	T_w=600 #[K] desired temperature of nozzle 


	## CONSTANTS OF SIM
	alpha = 0.07/8 # 0.07/8 : 1 ratio of alpha : beta gives very similar weights
	beta = 1
	design_alt = 9144
	truncate_ratio = 0.2# bounds on truncate < 0.1425

	CEA = CEA_constants(0) # not a functioning class as of now



	optimizer = aerospike_optimizer(r_e,T_w,alpha,beta,design_alt,truncate_ratio,chr_mesh_n=200,no_alt_range = 600,no_core=4)

	optimizer.multicore_animate_range(downstream_factor=1.2,chr_mesh_n=200,no_core=4)
	# contours = np.concatenate((optimizer.spike_opt.x,optimizer.spike_opt.y))

	# print(optimizer.cost_opt_contour_params(contours))

	#print(optimizer.cost_opt_design_params([design_alt,truncate_ratio]))
	#plt.legend()

	#optimizer.multicore_animate_to_its(812,4)


	## NEXT JOB, CHANGE TRUNCATION LENGTH TO 0.2, PLOT MACH NUMBER
# def cost_opt_design_params(params) : return cost_func_design_params(params,T_w,CEA,r_e,alpha,beta,chr_mesh_n=30,no_core=4)

# print(cost_opt_design_params([design_alt,truncate_ratio]))
# print(cost_opt_contour_params(np.concatenate((spike_opt.x,spike_opt.y))))

#res = scipy.optimize.minimize(optimizer.cost_opt_contour_params,contours)#,constraints = cons)

# with open('spike_opt.pkl','wb') as output:
# 	pickle.dump(spike_opt,output,pickle.HIGHEST_PROTOCOL)

# with open('spike_points.pkl','wb') as output:
# 	pickle.dump(res,output,pikcle.HIGHEST_PROTOCOL)

# with open('meshes.pkl','rb') as input:
# 	meshes = pickle.load(input)


## CONVERTING TO OPTIMIZABLE FUNCTION

# def min_design_alt(X):
# 	return X[0] - 3000

# def max_design_alt(X):
# 	return -X[0] + 12000

# def min_truncate(X):
# 	return X[1] - 0.2

# def max_truncate(X):
# 	return -X[1] + 1

# cons = [{'type':'ineq', 'fun':min_design_alt},{'type':'ineq', 'fun':max_design_alt},{'type':'ineq', 'fun':min_truncate},{'type':'ineq', 'fun':max_truncate}]

# minmizer_kwargs = {"constraints":cons}

# res = scipy.optimize.basinhopping(cost_lambda,[design_alt,truncate_ratio],minimizer_kwargs=minmizer_kwargs)
# print(res)