# Import Libraries/Modules
from ReadFile import Read
import numpy as np
import astropy.units as u

def ComponentMass(filename, part_type):
    """
    This function takes a galaxy simulation file and a particle type 
        and returns the total mass of all those particles in the
        galaxy simulation.

    Inputs: filename (string), name of galaxy simulation file
            part_type(float), number representing the particle
                type whose mass will be returned.
                Halo Type = 1.0
                Disk Type = 2.0
                Bulge Type = 3.0

    Outputs: mass (float), total mass of the specified particle type
                in the galaxy simulation file (1e12 MSun).

    """

    # Read in the data of the file
    _, _, data = Read(filename)

    # Find all the data of the specified particle type
    index = np.where(data["type"]==part_type)[0]
    component_data = data[index]
    
    # Calculate the combined mass of the particle type
    mass = np.sum(component_data["m"]) # Units of 1e10 MSun

    # Convert from units of 1e10 MSun to 1e12 MSun
    mass = mass*(1e-2) 

    # Return the mass rounded to three decimal places
    return np.round(mass, 3)
