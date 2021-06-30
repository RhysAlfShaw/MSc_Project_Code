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
        self.theta = x[0] #where on the orbit inradians
        self.M_star = x[1]
        self.M_planet = x[2]
        self.phi = x[3]
        self.a = x[4]
        self.e = x[5]
        self.Np = len(x[2])
        self.theta_1 = x[6]
        self.phi = x[7]
        self.sig = x[8]

    def calculate_params(self):
        self.c = self.a*self.e
        self.b = np.sqrt((self.a**2)- (self.c**2))
        self.R_vec = np.zeros((self.Np,3))
        self.V_vec = np.zeros((self.Np,3))
        for i in range(0,self.Np):
            self.R_vec[i,0] = self.a[i]*np.cos(self.theta[i])-self.c[i]
            self.R_vec[i,1] = self.b[i]*np.sin(self.theta[i])
            self.R_vec[i,2] = 0
            self.R_vec[i,:] = self.rotation(self.R_vec[i,:],self.theta_1,self.phi,self.sig)  
            
            Vmax = np.sqrt(G*(self.M_star+self.M_planet[i])*M_sol * ((2)/(np.linalg.norm(self.R_vec[i])*AU) - (1)/(self.a[i]*AU)))
            self.V_vec[i,0] = Vmax*np.sin(self.theta[i])
            self.V_vec[i,1] = Vmax*np.cos(self.theta[i])
            self.V_vec[i,2] = 0
            self.V_vec[i,:] = self.rotation(self.V_vec[i,:],self.theta_1,self.phi,self.sig)

    def ellipse(self,t,a,b,c):
        x = -c + a*np.cos(t)
        y = b*np.sin(t)
        z = np.zeros(len(t))
        X = np.array([x,y,z])
        return  X.T

    def rotation(self,X,theta_1,phi,sig):
        M_rot = np.array([[np.cos(theta_1)*np.cos(phi),np.cos(theta_1)*np.sin(phi), - np.sin(theta_1)],
                [- np.cos(sig)*np.sin(phi) + np.sin(sig)*np.sin(theta_1)*np.cos(phi),np.cos(sig)*np.cos(phi) + np.sin(sig)*np.sin(theta_1)*np.sin(phi),np.sin(sig)*np.cos(theta_1)],
                [np.sin(sig)*np.sin(theta_1) + np.cos(sig)*np.sin(theta_1)*np.cos(phi),-np.sin(sig)*np.cos(theta_1) + np.cos(sig)*np.sin(theta_1)*np.sin(phi),np.sin(sig)*np.cos(theta_1)]])
        return  np.matmul(M_rot,X)
    
    def graph(self,dimension):
        t = np.linspace(0,2*np.pi, 200)
        if dimension == '2D':
            plt.figure(figsize=[5*np.max(self.a),5*np.max(self.b)]) #figsize changes to appropriate deimensions.
            for i in range(0,self.Np):
                X = self.ellipse(t,self.a[i],self.b[i],self.c[i])
                for j in range(0,200):
                    X[j] = self.rotation(X[j],self.theta_1,self.phi,self.sig)
                plt.plot(X[:,0],X[:,1],linestyle='--')
                plt.scatter(self.R_vec[i,0],self.R_vec[i,1],label='Planet'+str(i+1)+'M ='+str(self.M_planet[i]))
            plt.scatter(0,0,label='Star')
            plt.xlabel('X')
            plt.ylabel('Y')
            plt.legend()
            plt.show()

        if dimension == '3D':

            fig = plt.figure()
            ax = fig.add_subplot(111,projection='3d')
            ax.scatter(0,0,0,label='star',s=40)
            for i in range(0,self.Np):
                X = self.ellipse(t,self.a[i],self.b[i],self.c[i])
                for j in range(0,200):
                    X[j] = self.rotation(X[j],self.theta_1,self.phi,self.sig)
                ax.plot(X[:,0],X[:,1],X[:,2],linestyle='--')
                ax.scatter(self.R_vec[i,0],self.R_vec[i,1],self.R_vec[i,2],label='Planet'+str(i+1)+'M ='+str(self.M_planet[i]))
            ax.set_box_aspect([1*np.max(self.a),1*np.max(self.b),1]) #figsize changes to appropriate dimentions
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            plt.legend()
            plt.show()


