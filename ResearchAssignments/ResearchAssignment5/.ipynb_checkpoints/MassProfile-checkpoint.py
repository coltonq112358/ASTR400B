# Import Libraries
from CenterOfMass import CenterOfMass
from ReadFile import Read
import numpy as np
import astropy.units as u
from astropy.constants import G

# Convert the gravitational constant to the correct units
G = G.to(u.kpc*u.km**2/u.s**2/u.Msun)

class MassProfile:
    def __init__(self, galaxy, snap):
        """
        Class to calculate the Mass Profile of a Galaxy at a specific snapshot.

            PARAMETERS
            ----------
            galaxy: 'str'
                Specifies which galaxy we are analyzing.
                e.g. 'MW', 'M31', or 'M33
            snap: 'int'
                Number which specifies the snapshot number of the galaxy to analyze
                e.g. 0, 1, etc.
        """
        # Store the name of the galaxy
        self.gname = galaxy
        
        # add a string of the filenumber to the value "000"
        ilbl = "000"+str(snap)
        # Remove all but the last 3 digits
        ilbl = ilbl[-3:]
        self.filename = "%s_"%(galaxy) + ilbl + '.txt'

        # read data in the given file using Read
        self.time, self.total, self.data = Read(self.filename)

        # store the mass and positions of the particles
        self.m = self.data['m']
        self.x = self.data['x'] * u.kpc
        self.y = self.data['y'] * u.kpc
        self.z = self.data['z'] * u.kpc

    def TotalMass(self, ptype):
        """
        Function to determine the total mass of a specific particly type
        within the galaxy.

        Inputs:
            ptype: Number indicating the particle type
                   Halo Type = 1.0
                   Disk Type = 2.0
                   Bulge Type = 3.0 
        
        Outputs:
            total_mass: Total mass of that particle type in MSun
                        (astropy quantity)
        """
        # Filter all particles by the specified particle type
        index = np.where(self.data['type']==ptype)
        # Calculate the total mass of that particle type
        total_mass = np.sum(self.m[index])
        # Return the total mass with the correct units
        return total_mass*1e10*u.Msun

    def MassEnclosed(self, ptype, radii):
        """
        Function which determines the mass enclosed within the galaxy for a list of radii

        Inputs:
            pytpe: number indicating the particle type
                   Halo Type = 1.0
                   Disk Type = 2.0
                   Bulge Type = 3.0
                   
            radii: numpy array of various radii
                   (kpc)

        Outputs:
            masses: numpy array of astropy quantities of the mass enclosed 
                    within the radii specified in the input
                    (MSun)
        """

        # Center Of Mass of the galaxy
        # Using the disk particles to calculate the center of mass
        COM = CenterOfMass(self.filename, 2)

        # Center of Mass Position
        COM_P = COM.COM_P(0.1)

        # Find all the particles of the specified type
        index = np.where(self.data['type']==ptype)

        # Find the masses and x, y, z, positions of all
        # the particles of the specified type
        pm = self.m[index]
        px = self.x[index] - COM_P[0]
        py = self.y[index] - COM_P[1]
        pz = self.z[index] - COM_P[2]
        pdist = np.sqrt(px**2 + py**2 + pz**2)

        # Initialize an array of the same size as radii
        # but whose values are all 0
        masses = np.zeros_like(radii)

        # Loop through the radii array and find the 
        # mass enclosed within that radius in the galaxy
        for i, r in enumerate(radii):
            dist_index = np.where(pdist < r*u.kpc)

            mass_enclosed = np.sum(pm[dist_index])

            masses[i] = mass_enclosed

        # Return the masses with the correct units
        return masses*1e10*u.Msun

    def MassEnclosedTotal(self, radii):
        """
        Function to calculate the total mass enclosed at specified radii

        inputs:
            radii: Numpy array of various radii
                   (kpc)

        Outputs:
            TotalMass: Numpy array of masses of the total mass
                       of the galaxy enclosed within a radius
                       (Array<Astropy.Quantities MSun>) 
        """
        # M33 Does not have Bulge mass so we only add up HaloMass
        # and DiskMass for M33
        if self.gname == "M33":
            # Calculate the Mass Components and add them together
            HaloMass = self.MassEnclosed(1, radii)
            DiskMass = self.MassEnclosed(2, radii)
            TotalMass = HaloMass + DiskMass
        else:
            # Calculate the Mass Components and add them together
            HaloMass = self.MassEnclosed(1, radii)
            DiskMass = self.MassEnclosed(2, radii)
            BulgeMass = self.MassEnclosed(3, radii)
            TotalMass = HaloMass + DiskMass + BulgeMass
        # Return the Total Enclosed Mass of the galaxy
        return TotalMass

    def HernquistMass(self, radius, a, MHalo):
        """
        Function to calcualte the mass enclosed within a certain radius
        assuming a Hernquist mass profile.

        Inputs:
            radius: float radius to calculate the enclosed mass
                    (kpc)
            a:      float Hernquist Scale
                    (kpc)
            MHalo: float Total dark matter mass of the galaxy
        """
        # Calculate the mass enclosed within the radius
        # using a Hernquist mass profile
        M = MHalo*(radius**2)/((a+radius)**2)
        # Return the mass
        return M

    def CircularVelocity(self, ptype, radii):      
        """
        Function to calculate the circular velocity of a particular
        particle type at a specific radius.

        Inputs:
            ptype: Number indicating the particle type
                   Halo Type = 1.0
                   Disk Type = 2.0
                   Bulge Type = 3.0 
            radii: Specific radius to calculate the circular orbital velocity of
                   (kpc)
        
        Outpus:
            speeds: Circular Orbital Velocities at the specified radii
                    (km/s)
        """
        # Calculate the mass enclosed within the radius of that particle type
        enclosed_masses = self.MassEnclosed(ptype, radii)
        # Assuming a circular orbit calculate the orbital velocity at that radius
        speeds = np.sqrt(G*enclosed_masses/(radii*u.kpc))
        # Return the speeds rounded to 2 decimal places
        return np.round(speeds, 2)

    def CircularVelocityTotal(self, radii):
        """
        Function to calculate the circular orbital velocity
        using the total mass of the galaxy enclosed within a radius

        Inputs:
            radii: numpy array of radii to calculate the circular orbitla velocity of
                   (kpc)
        
        Outputs:
            speeds: numpy array of circular velocity values rounded to two decimals
                    (km/s)
        """
        # Determine the total mass enclosed at the radii
        TotalMass = self.MassEnclosedTotal(radii)
        # Calculate the orbital velocity using the mass enclosed within the radii
        speeds = np.sqrt(G*TotalMass/(radii*u.kpc))
        # Return the speeds rounded to two decimals
        return np.round(speeds, 2)

    def HernquistVCirc(self, radius, a, MHalo):
        """
        Function to calculate the circular orbital velocity
        at the specific radius assuming a Hernquist Mass Profile 

        Inputs:
            radius: float radius to calculate orbital vleocity at
                    (kpc)
            a:      Hernquist scale
                    (kpc)
            MHalo: Total dark matter mass of the galaxy
                    (MSun)
        
        Outputs:
            HVCirc: Circulat orbital velocities rounded to 2 decimals
        """
        # Calculate the mass enclosed within the radius
        # using the Hernquist mass profile
        HMass = self.HernquistMass(radius, a, MHalo)
        # Calculate the orbital velocity using the enclosed mass
        HVCirc = np.sqrt(G*HMass/(radius*u.kpc))
        # Return the velocities rounded to 2 decimals
        return np.round(HVCirc, 2)     