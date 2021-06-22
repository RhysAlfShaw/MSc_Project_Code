import numpy as np
import matplotlib.pyplot as plt
N = 1000             #(number of stars)
Rmax = 1000          #(AU)
                     #Defalut units also include Mass in solar masses.
X = np.zeros(N)
Y = np.zeros(N)
Z = np.zeros(N)

for i in range(0,N):
    x = np.random.uniform(-Rmax,Rmax)
    y = np.random.uniform(-Rmax,Rmax)
    z = np.random.uniform(-Rmax,Rmax)
    if np.sqrt(x**2+y**2+z**2) < Rmax:
        X[i] = x    
        Y[i] = y
        Z[i] = z    
    
Vx = np.zeros(N)
Vy = np.zeros(N)
Vz = np.zeros(N)

for i in range(0,N):
    Vx[i] = np.random.uniform(0,1)
    Vy[i] = np.random.uniform(0,1)
    Vz[i] = np.random.uniform(0,1)

   

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(X,Y,Z)
plt.show()
