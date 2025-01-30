# Import libraries
import numpy as np
import astropy.units as u

# This function takes a data file and returns
# The time of the simulation, number of particles
# and the data for each of the particles
# Inputs: filename (name of the simulation file)
# Returns: time of simulation, count of particles
#          particle data
def Read(filename):
    file = open(filename, 'r')   # Open the file
    line1 = file.readline()      # Read the first line
    label, value = line1.split() # Split the first line
    time = float(value)*u.Myr    # Determine the time in Myr

    line2 = file.readline()      # Read the second line
    label, value = line2.split() # Split the second line
    count = int(value)           # Determine the number of particles

    file.close()                 # Close the file

    # Make a numpy array from the rest of the particle data
    data = np.genfromtxt(filename, dtype=None, names=True, skip_header=3)

    return time, count, data
