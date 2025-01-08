import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

data0 = np.loadtxt(r'./orbit_cdm.txt')
data1 = np.loadtxt(r'./orbit_sidm.txt')

#N=100 

x0 = data0[:,2]
y0 = data0[:,3]

x1 = data1[:,2]
y1 = data1[:,3]

fig = plt.figure()
#plt.contour(xi, yi, zi,50,linewidths=(1.5))
#plt.colorbar();
plt.xlabel("$x [Mpc]$")
plt.ylabel("$y [Mpc]$")
#plt.xlim([-0.4, 0.4])
#plt.ylim([-0.2, 0.2])

plt.plot(x1, y1, '--o',markersize = 2,color='blue')
plt.plot(x0, y0, '-o',markersize = 2,color='red')
#plt.xscale('log')
plt.legend(['SIDM 2','CDM'])

plt.show()
fig.savefig('demo_orbit.eps')
fig.savefig('demo_orbit.png')




