# add_material_dict.py

import re
import numpy as np
import pandas as pd
import CifFile
import mantid as mtd
from mantid.simpleapi import CreateSampleWorkspace, SetSample
from mantid.kernel import Material

# Define uc_volume() function
def uc_volume(a, b, c, alpha, gamma, beta):
    # Find unit cell volume in A^3
    unit_cell_volume = a * b * c * np.sqrt(1 - np.cos(alpha * np.pi / 180)**2 - np.cos(beta * np.pi / 180)**2 - np.cos(gamma * np.pi / 180)**2 + 2*np.cos(alpha * np.pi / 180) * np.cos(beta * np.pi / 180) * np.cos(gamma * np.pi / 180))
    return unit_cell_volume

# Define read_cif() function
def read_cif(filepath):
    # Read Crystallographic Information File
    cif = CifFile.ReadCif(filepath)

    # Define keys of interest
    keys = ["_chemical_formula_sum", "_cell_length_a", "_cell_length_b", "_cell_length_c", "_cell_angle_alpha", "_cell_angle_beta", "_cell_angle_gamma", "_cell_volume", "_cell_formula_units_Z"]

    # Initialize empty dictionary
    cif_dict = {}

    # Iterate over each block in CIF
    for database_code_PCD, block in cif.items():
        # Print Pearson's Crystal Data browser code
        # print(f"Block: {database_code_PCD}")
    
        # Initialize empty dictionary
        block_data = {}

        # Iterate through keys of interest
        for i in keys:
            # Check key exists in block
            if i in block:
                # Store key in dictionary
                block_data[i] = block[i]
            # If key does not exist in block, print message
            else:
                # Store None in dictionary
                block_data[i] = None
                print(f"{i}: Not found in {database_code_PCD}.")

        # Add block_data dictionary to information dictionary
        cif_dict[database_code_PCD] = block_data

    # Check for index troubleshooting
    # for database_code_PCD, block_data in cif_dict.items():
        # Print Pearson's Crystal Data browser code
        # print(f"Block: {database_code_PCD}")

    # Print value for each key of interest
    # for key, value in block_data.items():
        # print(f"{key}: {value}")
        # print()

    # Define remove_parentheses() function
    def remove_parentheses(i):
        # Replace text within parentheses (including parentheses) with empty string
        # \( matches an opening parenthesis
        # \) matches a closing parenthesis
        # [^)]* matches any character except a closing parenthesis
        return re.sub(r"\([^)]*\)", "", i) if i else None

    # Extract variables of interest
    a = remove_parentheses(cif_dict[database_code_PCD]["_cell_length_a"])
    b = remove_parentheses(cif_dict[database_code_PCD]["_cell_length_b"])
    c = remove_parentheses(cif_dict[database_code_PCD]["_cell_length_c"])
    alpha = remove_parentheses(cif_dict[database_code_PCD]["_cell_angle_alpha"])
    beta = remove_parentheses(cif_dict[database_code_PCD]["_cell_angle_beta"])
    gamma = remove_parentheses(cif_dict[database_code_PCD]["_cell_angle_gamma"])
    Z_param = cif_dict[database_code_PCD]["_cell_formula_units_Z"]
    formula = cif_dict[database_code_PCD]["_chemical_formula_sum"]
    cell_volume = cif_dict[database_code_PCD]["_cell_volume"]

    # Find sample number density
    # Number density of sample in number of atoms or formula units per cubic Angstrom
    sample_n_density = eval(Z_param) / eval(cell_volume)

    # Split string by spaces
    formula_split = formula.split() if formula else []

    # Initialize empty dictionary
    elements = {}

    # Define regular expression pattern for dictionary
    # Group 1: [A-Za-z]+ matches any uppercase (A-Z) or lowercase (a-z) letter or sequence of letters (i.e., element)
    # Group 2: (\([\w\d]+\)) matches any sequence denoted within parentheses (i.e., isotope like (Li7))
    # Group 3: (\d+(\.\d+)?) matches any digit (0-9), possibly with a decimal
    pattern = re.compile(r"([A-Za-z]+|\([\w\d]+\))(\d+(\.\d+)?)")

    # Iterate through each element; return both index (i) and string for each element in chemical formula
    for i, string in enumerate(formula_split):
        # Match string against regular expression pattern
        match = pattern.match(string)
    
        # Check object is not None
        if match:
            # Extract element/Group 1 match
            element_symbol = match.group(1)

            # Extract subscript/Group 2 match, if None (no subscript) then 1
            # None is unlikely when retrieving from "_chemical_formula_sum" in CIF
            subscript = match.group(2) if match.group(2) else "1"

            # Find how many numbers per formula units
            n_fu = eval(Z_param) * eval(subscript)

            # Create new entry in dictionary for each element
            elements[f"element_{i+1}"] = {"symbol": element_symbol, "subscript": subscript, "n_fu": n_fu}

    # Print dictionary for troubleshooting
    # for key, value in elements.items():
        # print(f"{key}: {value}")

    # Initialize empty array
    array = []

    # Iterate over each item in dictionary
    for key, value in elements.items():
        # Retrieve element symbol and numbers per formula units 
        symbol = value["symbol"]
        n_fu = value["n_fu"]

        # Concatenate and append to array
        concatenate = f"{symbol}{n_fu}"
        array.append(concatenate)

    # Print array for troubleshooting
    # print(array)

    # Concatenate array to string and replace spaces with "-"
    mantid_formula = "-".join([str(i) for i in array])
    # mantid_formula = mantid_formula.replace(" ", "-")

    # Initialize variable
    total_n = 0

    # Iterate over dictionary items
    for key, value in elements.items():
        # Extract numbers per formula units value and add to total sum
        n_fu = value["n_fu"]
        total_n += n_fu

    # Print total n
    # print("Total n:", total_n)
    
    # Return mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma
    return mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma, Z_param

# Define add_material_dict() function
def add_material_dict(cif_filepath, material_csv_filepath):
    # Extract mantid_formula, a, b, c, alpha, beta, gamma with read_cif() function
    mantid_formula, _, _, a, b, c, alpha, beta, gamma, Z_param = read_cif(cif_filepath)

    # Create empty data container/workspace
    ws = CreateSampleWorkspace()

    # Find unit cell volume in A^3
    unit_cell_volume = uc_volume(float(a), float(b), float(c), float(alpha), float(beta), float(gamma))

    # Add material to data container/workspace
    SetSample(ws, Material={"ChemicalFormula": mantid_formula,
                            "UnitCellVolume": float(unit_cell_volume),
                            "ZParameter": float(Z_param)})

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
        return re.sub(r"\d+", "", i)

    # Create new_row DataFrame 
    new_row = {"material": remove_subscript(mantid_formula),
               "absorb_xs": absorb_xs,
               "scatter_xs": scatter_xs,
               "scatter_length": scatter_length,
               "molecular_mass": molecular_mass,
               "unit_cell_volume": unit_cell_volume,
               "Z_param": Z_param,
    }

    # Insert new_row DataFrame to end of material_df DataFrame and reset index
    material_df.loc[len(material_df)] = new_row
    material_df = material_df.reset_index(drop = True)

    # Write material_df as csv to material_csv_filepath
    material_df.to_csv(material_csv_filepath, index = False)