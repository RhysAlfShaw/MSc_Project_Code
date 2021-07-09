from scipy.io import FortranFile as FortranFile
import numpy as np

class data:
    def __init__(self,filename):
        self.filename = filename


    def read(self):
        f = FortranFile('Simple_planet_test.dat', 'r')
        # number of bodies
        self.N = f.read_ints(np.int32)[0]
        print('N = ', self.N)
        # time 
        self.time = f.read_reals(np.double)[0]
        print('time', self.time)

        # unit system
        units = f.read_reals(np.double)
        self.umass = units[0]
        self.udist = units[1]
        self.utime = units[2]
        print('units: ', self.umass, self.udist, self.utime)
  
        self.x = f.read_reals(np.double)
  
        self.y = f.read_reals(np.double)
  
        self.z = f.read_reals(np.double)
  
        self.vx = f.read_reals(np.double)
  
        self.vy = f.read_reals(np.double)
  
        self.vz = f.read_reals(np.double)
 
        self.mass = f.read_reals(np.double)
        print('mass array', self.mass) 
        f.close()

        print(self.filename,' Data Read successfully!')
    
    def write(self,R_,V_,Nsys,t_start,Mass,Mass_unit=1.9884699E33,Length_unit=3.0856775E18,time_unit=31556952):
        #setting up the data to be written correctly.
        X = R_[:,0]
        Y = R_[:,1]
        Z = R_[:,2]
        Vx = V_[:,0]
        Vy = V_[:,1]
        Vz = V_[:,2]
        N = np.array([Nsys],dtype=np.int32)
        time = np.array([t_start],dtype=np.double)
        units = np.array([Mass_unit,Length_unit,time_unit],dtype=np.double)

        ## Following code writes the data in the required way.

        Data = FortranFile(self.filename,"w")
        Data.write_record(N)
        Data.write_record(time)
        Data.write_record(units)
        Data.write_record(X)
        Data.write_record(Y)
        Data.write_record(Z)
        Data.write_record(Vx)
        Data.write_record(Vy)
        Data.write_record(Vz)
        Data.write_record(Mass)
        Data.close()
        print(self.filename,' sucessfully written to directory!')