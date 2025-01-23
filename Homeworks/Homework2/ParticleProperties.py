import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, particle_type, particle_num):
    _, _, data = Read(filename)
    index = np.where(data['type'] == particle_type)[0]

    x = data['x'][index][particle_num]
    y = data['y'][index][particle_num]
    z = data['z'][index][particle_num]

    dist = np.sqrt(x**2 + y**2 + z**2)*u.kpc

    vx = data['vx'][index][particle_num]
    vy = data['vy'][index][particle_num]
    vz = data['vz'][index][particle_num]

    vel = np.sqrt(vx**2 + vy**2 + vz**2)*u.km/u.s

    mass = data['m'][index][particle_num]*u.solMass

    return dist, vel, mass

filename = "C:\\Users\\colto\\OneDrive\\Desktop\\ASTR400B\\Homeworks\\Homework2\\MW_000.txt"
dist, vel, mass = ParticleInfo(filename, 2.0, 99)
print(dist)
print(vel)
print(mass)