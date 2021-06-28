import numpy as np
import matplotlib.pyplot as plt

"""
Constants and Unit conversions
"""

M_sol = 1.99E30
G = 6.67E-11
AU = 1.49E11



class gen_cluster:
    """
    This class will create a cluster, where options exist for different star distributions
    different IMFs and for other parameters such as N, radius of clusters sphere,
    and the virial ratio q.
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
        if IMF == 'KROUPA':
            self.Mass = self.cust_dis(0,20,self.Kroupa_IMF)
        else:
            print("INPUT ERROR: Please specify 'constant' or....")

    def cust_dis(self,x0,x1,imf,nControl=10**6):
        sample = []
        nLoop  = 0
        while len(sample)<self.N and nLoop<nControl:
            x = np.random.uniform(x0,x1)     
            prop = imf(x)
            assert prop>=0
            if np.random.uniform(0,1) <= prop: #26.67 is max probability 
                sample += [x]
            nLoop+= 1
        return np.array(sample)
    
    def Kroupa_IMF(self,m):
        alpha = 0
        if m > 0 and m < 0.08:
            alpha = 0.3
        if m > 0.08 and m < 0.5:
            alpha = 1.3
        if m > 0.5:
            alpha = 2.3
        return m**(-alpha)


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

            plt.figure(figsize=[3*np.max(self.a),3*np.max(self.b)]) #figsize changes to appropriate deimensions.
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


        




class measure:

    """
    allows for measurements of a system. Giving only R_vec and V_vec. This will be important
    when measuring the properties of otbits over time.
    
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