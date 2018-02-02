import numpy as np
import matplotlib.pyplot as plt
import multiprocessing as mp
import pickle
from plug_nozzle_angelino import plug_nozzle
from MOC import chr_mesh

r_e = 0.067/2 #0.034 # likely too large
expansion_ratio = 8.67 #8.1273
A_t = r_e**2*np.pi/expansion_ratio # max expansion (r_b = 0, r_e**2 >= A_t*expansion_ratio/np.pi)
gamma = 1.2381 #np.mean([1.2534,1.2852])
T_c = 2833.63
p_c = 34.474*10**5
rho_c = 3.3826
R = (1-1/gamma)*1.8292*1000#1.8292 1.6196
a_c = np.sqrt(gamma*(1-1/gamma)*200.07*T_c) 
n = 1000

def thrust_curve(spike,gamma,altitude,num_wave,downstream_factor,send_end):

	thrust = [] #np.zeros(num_wave.shape)

	for i in range(len(num_wave)):
		#print(num_wave[i])
		MOC_mesh = chr_mesh(spike,gamma,altitude,num_wave[i],downstream_factor=downstream_factor,plot_chr=0)
		
		thrust.append(MOC_mesh)

		##
		#thrust[i] = MOC_mesh.compute_thrust('cubic',100)

	send_end.send(MOC_mesh)


def multicore_thrust_compute(spike,gamma,altitude,num_wave,downstream_factor,no_core):
	proc_list = []
	pipe_list =[]
	np.random.shuffle(num_wave)
	num_wave_split = np.split(num_wave,no_core)
	
	for i in range(no_core):
		recv_end, send_end = mp.Pipe(False)
		args = (spike,gamma,altitude,num_wave_split[i],downstream_factor,send_end)

		proc = mp.Process(target=thrust_curve, args = args)
		proc_list.append(proc)
		pipe_list.append(recv_end)
		proc.start()

	for proc in proc_list:
		proc.join()

	thrust = [x.recv() for x in pipe_list]

	print(thrust)

	thrust = np.concatenate(thrust_range)

	# for thread in threads:
	# 	thread.map

	return thrust

no_core = 8

num_wave = np.arange(25*no_core); #num_wave = num_wave[29:]

num_wave = num_wave[4*no_core:]

altitude = 6000
# design plug nozzle


spike = plug_nozzle(expansion_ratio,A_t,r_e,gamma,T_c,p_c,a_c,rho_c,n,truncate_ratio = 0.2)

print("Computing meshes with " + str(no_core) + " cores")
# meshes = multicore_thrust_compute(spike,gamma,altitude,num_wave,1.2,no_core)

mesh = chr_mesh(spike,gamma,6000,50,downstream_factor=1.2,plot_chr=0)
print(mesh.compute_thrust('nearest',30))
contour_line, = plt.plot(spike.x,spike.y,'b', label = 'Nozzle Contour')
centre_line, = plt.plot(spike.x,np.zeros(spike.x.shape),'k--',label = 'Centre Line')
lip_coord, = plt.plot(spike.lip_x,spike.lip_y,'rx', label = 'Lip Location')



plt.legend(handles = [contour_line,centre_line,lip_coord])
plt.xlabel('(m)')
plt.ylabel('(m)')
print(mesh.PR)
#plt.title('Sea Level (PR = ' + str(np.round(mesh.PR,0)) + ') Characteristic Mesh')
plt.axis('equal')
plt.show()

#mesh.compute_thrust('nearest',30)
# print("Pickling result")
# with open('meshes.pkl','wb') as output:
# 	pickle.dump(meshes,output,pickle.HIGHEST_PROTOCOL)

# with open('meshes.pkl','rb') as input:
# 	meshes = pickle.load(input)


# plt.plot(num_wave,thrust,'bo')


#plt.axis('equal')
