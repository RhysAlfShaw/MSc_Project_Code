import numpy as np
import matplotlib.pyplot as plt

"""
Constants and Unit conversions
"""
M_planet = 9.55E-4 #M_sol Jupiter like mass
M_sol = 1.99E30
G = 6.67E-11
AU = 1.49E11

def ellipse(a,b,t):
    u = -c     #x-position of the center
    x = u + a*np.cos(t)
    y = b*np.sin(t)
    return x,y

class gen_cluster:
    """
    This class will create 
    """
    def __init__(self,x):
         self.i = x**2

class get_planet:
    """
    This class will create a planetary systems based on input parameters on a clsuters system.
    """
    def __init__(self):
        self.i = 0



class measure:
    """
    allows for measurements of a system.
    """
    def __init__(self,data):
        self.M_star = data[0]
        self.M_planet = data[1]
        self.R_vec = data[2]
        self.V_vec = data[3]
    
    def calculate_a(self):
        mu = self.M_star*M_sol*M_sol*self.M_planet/(self.M_planet*M_sol+self.M_star*M_sol)
        Eb = abs(0.5*mu*np.linalg.norm(self.V_vec)**2  - (G*self.M_star*M_sol*self.M_planet*M_sol)/(np.linalg.norm(self.R_vec)*AU))
        a_calc = (G * self.M_star*M_sol * self.M_planet*M_sol) / (2*Eb*AU)
        return a_calc #units of AU
    
    def calculate_e(self,a):
        e = np.sqrt((1-(np.linalg.norm(self.R_vec))/(abs(a)) )**2 + (np.dot(self.R_vec,self.V_vec)**2)/(abs(a)*G*(self.M_planet+self.M_star)*M_sol))
        return e
    
    def graph(self,dimension,a,e):
        c = a*e #from the definition of eccentricity.
        b = np.sqrt((a**2) - (c**2))
        t = np.linspace(0,2*np.pi, 200)
        x = -c + a*np.cos(t)
        y = b*np.sin(t)
        if dimension == '2D':
            plt.figure(figsize=[5*a,5*b]) #figsize changes to appropriate deimensions.
            plt.plot(x,y,linestyle='--')
            plt.scatter(0,0,label='Star')
            plt.scatter(self.R_vec[0],self.R_vec[1],label='Planet')
            plt.legend()
            plt.show()
        if dimension == '3D':
            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(0,0,0,label='star',s=40)
            ax.plot(x,y,0,linestyle='--')
            ax.scatter(self.R_vec[0],self.R_vec[1],self.R_vec[2],label='Planet')
            ax.set_box_aspect([1*a,1*b,1]) #figsize changes to appropriate dimentions
            plt.legend()
            plt.show()
        else:
            print("INPUT ERROR: Please specify '2D' or '3D' only!")