import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate

data0 = np.loadtxt(r'./outbind.txt')
#data1 = np.loadtxt(r'./outunbind.txt')

N=10 

x = data0[:,0]
y = data0[:,1]
z = data0[:,2]

#x1 = data1[:,0]
#y1 = data1[:,1]
#z1 = data1[:,2]

def rate(x,y):
    return z

xi = np.linspace(x.min(), x.max(), N)
yi = np.linspace(y.min(), y.max(), N)
zi = scipy.interpolate.griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')


fig = plt.figure()
plt.contour(xi, yi, zi,50,linewidths=(1.5))
plt.colorbar();
plt.xlabel("$x [Mpc]$")
plt.ylabel("$y [Mpc]$")
#plt.xlim([-0.4, 0.4])
#plt.ylim([-0.2, 0.2])

#plt.plot(x1, y1, 'o',markersize = 0.02,color='blue')
#plt.plot(x0, y0, 'o',markersize = 0.01,color='red')
plt.show()
#fig.savefig('demo_snap011_100m.eps')
#fig.savefig('demo_snap011_100m.png')




