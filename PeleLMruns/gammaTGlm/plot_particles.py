import os
import glob
import numpy as np
from matplotlib import pyplot as plt

def read_particle_ascii_data(filename):
    with open(filename,'r') as f:
        n_particles = int(f.readline())
        for i in range(4): f.readline()
        # gotta make this general at some point
        data = np.zeros((n_particles,10))
        for n,line in enumerate(f):
            data[n,:] = [float(i) for i in line.split()]

    return data


files = sorted(glob.glob('spray*.p3d'))
for file in files:
    data = read_particle_ascii_data(file)
    plt.scatter(data[:,0],data[:,1])
    plt.show()
    plt.close()
