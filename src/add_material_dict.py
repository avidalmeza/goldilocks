# add_material_dict.py

import re
import numpy as np
import pandas as pd
import mantid as mtd
from mantid.simpleapi import CreateSampleWorkspace, SetSample
from mantid.kernel import Material
from read_cif import read_cif

# Define uc_volume() function
def uc_volume(a, b, c, alpha, gamma, beta):
    # Find unit cell volume in A^3
    unit_cell_volume = a * b * c * np.sqrt(1-np.cos(alpha*np.pi/180)**2 - np.cos(beta*np.pi/180)**2 - np.cos(gamma*np.pi/180)**2 + 2*np.cos(alpha*np.pi/180) * np.cos(beta*np.pi/180) * np.cos(gamma*np.pi/180))
    return unit_cell_volume

# Define add_material_dict() function
def add_material_dict(cif_filepath, material_csv_filepath):
    # Extract mantid_formula, a, b, c, alpha, beta, gamma with read_cif() function
    mantid_formula, _, _, a, b, c, alpha, beta, gamma, Z_param = read_cif(cif_filepath)

    # Create empty data container/workspace
    ws = CreateSampleWorkspace()

    # Find unit cell volume in A^3
    unit_cell_volume = uc_volume(float(a), float(b), float(c), float(alpha), float(beta), float(gamma))

    # Add material to data container/workspace
    SetSample(ws, Material={'ChemicalFormula': mantid_formula,
                            'UnitCellVolume': float(unit_cell_volume),
                            'ZParameter': float(Z_param)})

    # Obtain sample object from workspace
    sample = ws.sample()

    # Retrieve absorption cross-section
    # Define absorption cross section per formula unit in bn/fu
    absorb_xs = float(sample.getMaterial().absorbXSection())

    # Retrieve total scattering cross-section
    # Define scattering cross section per formula unit in bn/fu
    scatter_xs = float(sample.getMaterial().totalScatterXSection())

    # Retrieve coherent scattering length in fm
    scatter_length = float(sample.getMaterial().cohScatterLength())

    # Retrieve relative molecular mass
    # Define molecular mass in g/mol/fu
    molecular_mass = float(sample.getMaterial().relativeMolecularMass())

    # Read sample can material csv
    material_df = pd.read_csv(material_csv_filepath)

    # Define remove_subscript() function
    def remove_subscript(i):
        # Replace all digits in string with empty string
        return re.sub(r'\d+', '', i)

    # Create `new_row` DataFrame 
    new_row = {'material': remove_subscript(mantid_formula),
               'absorb_xs': absorb_xs,
               'scatter_xs': scatter_xs,
               'scatter_length': scatter_length,
               'molecular_mass': molecular_mass,
               'unit_cell_volume': unit_cell_volume,
               'Z_param': Z_param,
    }

    # Insert `new_row` DataFrame to end of `material_df` DataFrame and reset index
    material_df.loc[len(material_df)] = new_row
    material_df = material_df.reset_index(drop = True)

    # Write `material_df` as csv to material_csv_filepath
    material_df.to_csv(material_csv_filepath, index = False)