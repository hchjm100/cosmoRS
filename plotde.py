import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

data0 = np.loadtxt(r'./outdebug.txt')
data1 = np.loadtxt(r'./center.txt')

#N=100 

x0 = data0[:,0]
y0 = data0[:,1]
z0 = data0[:,2]

x1 = data1[:,0]
y1 = data1[:,1]

fig = plt.figure()
plt.xlabel("$x [Mpc]$")
plt.ylabel("$y [Mpc]$")
plt.xlim([-0.4, 0.4])
plt.ylim([-0.2, 0.2])

plt.plot(x0, y0, 'o',markersize = 0.08,color='red')
plt.plot(x1, y1, 'o',markersize = 3,color='blue')
plt.show()




