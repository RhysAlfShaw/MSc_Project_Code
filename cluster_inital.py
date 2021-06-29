import numpy as np
import matplotlib.pyplot as plt

import clustersimpy as cs

N = 1000              #(number of stars)
Rmax = 50000
q = 0.5
D = 1
x = [Rmax,N,q]
R = cs.gen_cluster(x)
R.Gen_positions('Fractal',D=D)
R.Gen_mass('KROUPA')
R.Gen_velocities()
R.graph('3D')
R.graph('2D')

