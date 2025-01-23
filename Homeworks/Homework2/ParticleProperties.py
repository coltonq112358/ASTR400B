import numpy as np
import astropy.units as u
from ReadFile import Read

def ParticleInfo(filename, particle_type, particle_num):
    _, _, data = Read(filename)

    particle_num = int(particle_num)

    x = data['x'][particle_num]
    y = data['y'][particle_num]
    z = data['z'][particle_num]

    dist = np.sqrt(x**2 + y**2 + z**2)*u.kpc

    vx = data['vx'][particle_num]
    vy = data['vy'][particle_num]
    vz = data['vz'][particle_num]

    vel = np.sqrt(vx**2 + vy**2 + z**2)*u.km/u.s

    mass = data['m'][particle_num]*u.solMass

    return dist, vel, mass

name = "MW_000.txt"
dist, vel, mass = ParticleInfo(name, 1.0, 1)
print(dist, vel, mass)



