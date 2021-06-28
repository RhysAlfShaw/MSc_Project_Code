import numpy as np
import matplotlib.pyplot as plt
import genplanet_test as gp
"""
N = 1000             #(number of stars)
Rmax = 1000          #(AU)
                     #Defalut units also include Mass in solar masses.
X = np.zeros(N)
Y = np.zeros(N)
Z = np.zeros(N)
i = 0
while i < N:
    x = np.random.uniform(-Rmax,Rmax)
    y = np.random.uniform(-Rmax,Rmax)
    z = np.random.uniform(-Rmax,Rmax)
    if np.sqrt(x**2+y**2+z**2) < Rmax:
        X[i] = x    
        Y[i] = y
        Z[i] = z
        i = i + 1    
    
Vx = np.zeros(N)
Vy = np.zeros(N)
Vz = np.zeros(N)

for i in range(0,N):
    Vx[i] = np.random.uniform(0,1)
    Vy[i] = np.random.uniform(0,1)
    Vz[i] = np.random.uniform(0,1)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.set_box_aspect([1,1,1])
scat = ax.scatter(X,Y,Z)
plt.show()
"""
theta = np.array([0, 0, 0, 0])
M_star = 0.2
M_planet = np.array([9.55E-4,9.55E-4,9.55E-4,9.55E-4])
phi = np.array([0,0,0,0])
a = np.array([1,3,5,6]) #semi_major axis of the orbits in AU
e = np.array([0,0,0,0]) # all circular orbits!

x = [theta,M_star,M_planet,phi,a,e]
P = gp.get_planet(x)
P.calculate_params()
P.graph('3D')
