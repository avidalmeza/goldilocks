import os
from goldilocks import xs_calculator

# Set sample filepath
# E.g., src/cif/Ba2Co1La2O12Te2.cif
thisDir = os.getcwd()
sample = os.path.join(thisDir, 'src', 'cif', 'Ba2Co1La2O12Te2.cif')

# Execute calculator
xs_calculator(x = sample, neutron_energy = 100, pack_fraction = 0.50)