import numpy as np
import matplotlib.pyplot as plt

"""
Constants and Unit conversions
"""

M_sol = 1.99E30
G = 6.67E-11
AU = 1.49E11       


class get_planet:
    """
    This class will create a planetary systems based on input parameters on a clsuters system.
    Takes the inputs of theta (to pick the starting position of the orbit), M_star, M_planet, Phi and eccentricity and semi major axis a.
    Phi should all for the inclination of the orbit to be changed. 
    """
    def __init__(self,x): #x must be a list of each of these parameters! 
        self.theta = x[0]
        self.M_star = x[1]
        self.M_planet = x[2]
        self.phi = x[3]
        self.a = x[4]
        self.e = x[5]
        self.Np = len(x[2])

    def calculate_params(self):
        self.c = self.a*self.e
        self.b = np.sqrt((self.a**2)- (self.c**2))
        self.R_vec = np.zeros((self.Np,3))
        self.V_vec = np.zeros((self.Np,3))
        for i in range(0,self.Np):
            self.R_vec[i,0] = self.a[i]*np.cos(self.theta[i])-self.c[i]
            self.R_vec[i,1] = self.b[i]*np.sin(self.theta[i])
            self.R_vec[i,2] = 0
            Vmax = np.sqrt(G*(self.M_star+self.M_planet[i])*M_sol * ((2)/(np.linalg.norm(self.R_vec[i])*AU) - (1)/(self.a[i]*AU)))
            self.V_vec[i,0] = Vmax*np.sin(self.theta[i])*np.cos(self.phi[i])
            self.V_vec[i,1] = Vmax*np.cos(self.theta[i])*np.cos(self.phi[i])
            self.V_vec[i,2] = Vmax*np.sin(self.phi[i])
            
    
    def graph(self,dimension):
        t = np.linspace(0,2*np.pi, 200)
        x = np.zeros((self.Np,200))
        y = np.zeros((self.Np,200))
        for i in range(0,self.Np):
            x[i,:] = -self.c[i] + self.a[i]*np.cos(t)
            y[i,:] = self.b[i]*np.sin(t)
        
        if dimension == '2D':

            plt.figure(figsize=[5*np.max(self.a),5*np.max(self.b)]) #figsize changes to appropriate deimensions.
            for i in range(0,self.Np):
                plt.plot(x[i],y[i],linestyle='--')
                plt.scatter(self.R_vec[i,0],self.R_vec[i,1],label='Planet'+str(i+1)+'M ='+str(self.M_planet[i]))
            plt.scatter(0,0,label='Star')
            plt.legend()
            plt.show()

        if dimension == '3D':

            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(0,0,0,label='star',s=40)
            for i in range(0,self.Np):
                ax.plot(x[i],y[i],0,linestyle='--')
                ax.scatter(self.R_vec[i,0],self.R_vec[i,1],self.R_vec[i,2],label='Planet'+str(i+1)+'M ='+str(self.M_planet[i]))
            ax.set_box_aspect([1*np.max(self.a),1*np.max(self.b),1]) #figsize changes to appropriate dimentions
            plt.legend()
            plt.show()


