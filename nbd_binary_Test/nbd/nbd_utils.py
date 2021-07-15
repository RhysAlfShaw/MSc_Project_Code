#!/usr/bin/python
import os
import numpy as np
from scipy.io import FortranFile as FortranFile

#
# define a class to hold all our nbd data. The class (like a struct) is
# what we return from the function
class nbd_data_struct:
    pass

def nbd_read(filename):
  #
  # empty class
  a = nbd_data_struct()  
  #
  # open file
  print('Opening NBD snapshot ', filename)
  f = FortranFile(filename, 'r')
  #
  # number of bodies
  a.N = f.read_ints(np.int32)[0]
  print('N = ', a.N)
  #
  # time 
  a.time = f.read_reals(np.double)[0]
  print('time', a.time)
  #
  # unit system
  units = f.read_reals(np.double)
  a.umass = units[0]
  a.udist = units[1]
  a.utime = units[2]
  print('units', a.umass, a.udist, a.utime)
  #
  # x pos
  a.x = f.read_reals(np.double)
  #
  # y pos
  a.y = f.read_reals(np.double)
  #
  # z pos
  a.z = f.read_reals(np.double)
  #
  # vx pos
  a.vx = f.read_reals(np.double)
  # 
  # vy pos
  a.vy = f.read_reals(np.double)
  #
  # vz pos
  a.vz = f.read_reals(np.double)
  #
  # mass
  a.mass = f.read_reals(np.double)
  print('mass array', a.mass) 
  #
  # close file!
  f.close()
  #
  # we simply return the struct a.--- that holds all the data! 
  return a  
