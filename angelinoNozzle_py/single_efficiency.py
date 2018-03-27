import numpy as np 
import matplotlib.pyplot as plt 
from scipy import interpolate

PRs = np.linspace(8.9,200,100)

s_pr = np.array([8.9,10,20,30,100,200])
s_eff = np.array([0.89,0.88,0.82,0.85,0.94,0.97])



b_pr = np.array([8.9,30,200,400])
b_eff = np.array([0.18,0.86,1,0.92])

average_b = np.trapz(b_eff[:-1],b_pr[:-1])/(b_pr[-2]-b_pr[0])
average_s = np.trapz(s_eff,s_pr)/(s_pr[-1]-s_pr[0])

print("Aerospike avg. eff: " + str(average_s))
print("Bell avg. eff: " + str(average_b))

plt.plot(s_pr,s_eff,'b-o',label = '20% Truncated Nozzle')
plt.plot(b_pr[:-1],b_eff[:-1],'-o',color='C1',label='Bell Nozzle')

plt.legend()

plt.title('Efficiency Currve for Bell and Truncated Aerospike')
plt.xlabel('PR')
plt.ylabel('Efficiency')
plt.show()