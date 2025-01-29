# Import libraries
import numpy as np
import astropy.units as u
from ReadFile import Read

# This function finds properties of a 
# particular particle
# Inputs: filename (name of simulation file),
#         particle_type (number specifying type of particle),
#         particle_num (number specifying which particle)
# Returns: dist (Magnitude of distance in kpc)
#          vel (Magnitude of velocity in km/s)
#          mass (Mass of the particle in solar Mass)
def ParticleInfo(filename, particle_type, particle_num):
    _, _, data = Read(filename)                        # Read the data from the file
    index = np.where(data['type'] == particle_type)[0] # Find the indices of the particle type

    # Determine the x, y, and z position of the particle
    x = data['x'][index][particle_num]
    y = data['y'][index][particle_num]
    z = data['z'][index][particle_num]

    # Determine the distance of the particle
    # Round to 3 decimal places
    dist = np.around(np.sqrt(x**2 + y**2 + z**2)*u.kpc, 3)
    
    # Determine the velocity of the particle in x, y, and z
    vx = data['vx'][index][particle_num]
    vy = data['vy'][index][particle_num]
    vz = data['vz'][index][particle_num]
    
    # Determine the magnitude of the velocity of the particle
    # round to 3 decimal places
    vel = np.around(np.sqrt(vx**2 + vy**2 + vz**2)*u.km/u.s, 3)

    # Determine the mass of the particle
    mass = data['m'][index][particle_num]*1e10*u.solMass

    return dist, vel, mass
