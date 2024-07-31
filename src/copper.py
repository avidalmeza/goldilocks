import numpy as np
import mantid as mtd
from mantid.simpleapi import CreateSampleWorkspace, SetSample
from mantid.kernel import Material

# Define unit_cell_volume() function
def uc_volume(a, b, c, alpha, gamma, beta):
    # Find unit cell volume in A^3
    unit_cell_volume = a * b * c * np.sqrt(1-np.cos(alpha * np.pi/180)**2 - np.cos(beta * np.pi/180)**2 - np.cos(gamma * np.pi/180)**2 + 2*np.cos(alpha * np.pi/180) * np.cos(beta * np.pi/180) * np.cos(gamma * np.pi/180))
    return unit_cell_volume

# Set a, b, c, alpha, beta, gamma, Z for common alumnium 
a = 3.6199
b = 3.6199
c = 3.6199
alpha = 90.0
beta = 90.0
gamma = 90.0
Z_param = 4.0

# Find unit cell volume in A^3
unit_cell_volume = uc_volume(a, b, c, alpha, beta, gamma)
print(unit_cell_volume)

# Create empty data container/workspace
ws = CreateSampleWorkspace()

# Add material to data container/workspace
SetSample(ws, Material = {'ChemicalFormula': 'Cu',
                          'UnitCellVolume': float(unit_cell_volume),
                          'ZParameter': Z_param})

# Obtain sample object from workspace
sample = ws.sample()

# Retrieve absorption cross-section
# Define absorption cross section per formula unit in bn/fu
Cu_absorb_xs = float(sample.getMaterial().absorbXSection())
print(f'Absorption cross-section: {Cu_absorb_xs}')

# Retrieve total scattering cross-section
# Define scattering cross section per formula unit in bn/fu
Cu_scatter_xs = float(sample.getMaterial().totalScatterXSection())
print(f'Total scattering cross-section: {Cu_scatter_xs}')

# Retrieve coherent scattering length
Cu_scatter_length = float(sample.getMaterial().cohScatterLength())
print(f'Coherent scattering length: {Cu_scatter_length}')

# Retrieve relative molecular mass
# Define molecular mass in g/mol/fu
Cu_molecular_mass = float(sample.getMaterial().relativeMolecularMass())
print(f'Relative molecular mass: {Cu_molecular_mass}')