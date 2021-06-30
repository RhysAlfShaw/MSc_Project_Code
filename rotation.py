import numpy as np 
import matplotlib.pyplot as plt 



# system parameters
G = 6.67E-11
M_star = 0.2       #M_sol
M_planet = 9.55E-4 #M_sol Jupiter like mass
M_sol = 1.99E30
AU = 1.49E11 
a = 1                #AU Semi-major axis
r = 0.8              #AU min distance between star and planet
c = a - r
b = np.sqrt((a**2) - (c**2))

def ellipse(a,b,t):
    u = -c     #x-position of the center
    x = u + a*np.cos(t)
    y = b*np.sin(t)
    z = np.zeros(len(t))
    X = np.array([x,y,z])
    return  X.T

t = np.linspace(0,2*np.pi, 200)
X = ellipse(a,b,t)
print(X[0])

theta = 0
R_vec = np.array([a*np.cos(theta)-c,b*np.sin(theta),0]) #-c important to account for shift in origin
print(R_vec)

Vmax = np.sqrt(G*(M_star+M_planet)*M_sol * ((2)/(np.linalg.norm(R_vec)*AU) - (1)/(a*AU)))
V_vec = np.array([Vmax*np.sin(theta),Vmax*np.cos(theta),0])

theta_1 = 0 #rotation around x axis
phi = 0 #rotation around z axis
sig = 0.2 #rotation around y axis
M_rot = np.array([[np.cos(theta_1)*np.cos(phi),np.cos(theta_1)*np.sin(phi), - np.sin(theta_1)],
                [- np.cos(sig)*np.sin(phi) + np.sin(sig)*np.sin(theta_1)*np.cos(phi),np.cos(sig)*np.cos(phi) + np.sin(sig)*np.sin(theta_1)*np.sin(phi),np.sin(sig)*np.cos(theta_1)],
                [np.sin(sig)*np.sin(theta_1) + np.cos(sig)*np.sin(theta_1)*np.cos(phi),-np.sin(sig)*np.cos(theta_1) + np.cos(sig)*np.sin(theta_1)*np.sin(phi),np.sin(sig)*np.cos(theta_1)]])

R_rot = np.matmul(M_rot,R_vec)
print(R_rot)
for i in range(0,200):
    X[i] = np.matmul(M_rot,X[i])



fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
ax.scatter(0,0,0,label='star',s=40)
ax.plot(X[:,0],X[:,1],X[:,2],linestyle='--')
ax.scatter(R_rot[0],R_rot[1],R_rot[2],label='Planet')
ax.set_box_aspect([1*a,1*b,1]) #figsize changes to appropriate dimentions
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
plt.legend()
plt.show()
