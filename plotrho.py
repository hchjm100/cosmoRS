import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

data0 = np.loadtxt(r'./orbit_cdm.txt')
data1 = np.loadtxt(r'./orbit_sidm.txt')


#N=100 

x0 = data0[:,0]
y0 = data0[:,1]
x1 = data1[:,0]
y1 = data1[:,1]

fig = plt.figure()
plt.xlabel(r"$t [Gyr]$")
plt.ylabel(r"$\rho [M_{\odot}/pc^3]$")
plt.xlim([0, 10])
#plt.ylim([-0.2, 0.2])

plt.plot(x1, y1, '--o',markersize = 2,color='blue')
plt.plot(x0, y0, '-o',markersize = 2,color='red')
#plt.xscale('log')
plt.legend(['SIDM 2','CDM'])
plt.show()

fig.savefig('demo_rho_t_df2.eps')
fig.savefig('demo_rho_t_df2.png')

