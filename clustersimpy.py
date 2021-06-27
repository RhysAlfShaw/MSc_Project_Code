import numpy as np
import matplotlib.pyplot as plt

"""
Constants and Unit conversions
"""
M_planet = 9.55E-4 #M_sol Jupiter like mass
M_sol = 1.99E30
G = 6.67E-11
AU = 1.49E11

class gen_cluster:
    """
    This class will create 
    """
    def __init__(self,x):
        self.Rmax = x[0]
        self.N = x[1]
        self.q = x[2]
    
    def Gen_positions(self,distribution):
        X = np.zeros(self.N)
        Y = np.zeros(self.N)
        Z = np.zeros(self.N)
        if distribution == 'Uniform':
            i = 0
            while i < self.N:
                x = np.random.uniform(-self.Rmax,self.Rmax)
                y = np.random.uniform(-self.Rmax,self.Rmax)
                z = np.random.uniform(-self.Rmax,self.Rmax)
                if np.sqrt(x**2+y**2+z**2) < self.Rmax :
                    X[i] = x    
                    Y[i] = y
                    Z[i] = z
                    i = i + 1
            self.X = X
            self.Y = Y
            self.Z = Z

        else:
            print("INPUT ERROR: Please specify 'Uniform'or ....")
    
    def Gen_mass(self,IMF):
        if IMF == 'constant':
            Mass = np.zeros(self.N)

            for i in range(0,self.N):
                Mass[i] = 0.2
            self.Mass = Mass
        if IMF == 'OTHER':
            #OTHER IMF CODE HERE!
            i = 0 # place holder code!
        else:
            print("INPUT ERROR: Please specify 'constant' or....")

    def Gen_velocities(self):
        Vx = np.zeros(self.N)
        Vy = np.zeros(self.N)
        Vz = np.zeros(self.N)

        for i in range(0,self.N): # Velocities between 0 and 1.
            Vx[i] = np.random.uniform(0,1)
            Vy[i] = np.random.uniform(0,1)
            Vz[i] = np.random.uniform(0,1)
        PE_tot = 0
        #try:
        for k in range(1,self.N): #note since 0 is the begining of the index.
             for j in range(0,k-1):
                R_res = np.sqrt((self.X[k]-self.X[j])**2 + (self.Y[k]-self.Y[j])**2 + (self.Z[k]-self.Z[j])**2)
                PE_tot += (G * self.Mass[j]*M_sol**2*self.Mass[k])/(R_res*AU) 
        
        #else:
        #    print("MASS OR POSITION ERROR!!, Please run Gen_position & Gen_mass first to prevent this error!!")
        
        print('Potential Energy',PE_tot)
        # Calculating Kinetic Energy!
        Ek_tot = 0
        for i in range(0,self.N):
            Ek_tot += 0.5 * self.Mass[i] *M_sol* (np.sqrt(Vx[i]**2 + Vy[j]**2 + Vz[i]**2))**2
        
        a = np.sqrt((self.q*PE_tot)/(Ek_tot))
        self.Vx = Vx*a
        self.Vy = Vy*a
        self.Vz = Vz*a

    def graph(self):
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        ax.scatter(self.X,self.Y,self.Z)
        ax.set_box_aspect([1,1,1])
        plt.show()

        
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