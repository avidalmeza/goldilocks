# goldilocks.py

import os
import re
import numpy as np
import pandas as pd
import CifFile
import mantid as mtd
from mantid.simpleapi import CreateSampleWorkspace, SetSample
from mantid.kernel import Material

# Obtain current working directory filepath
this_dir = os.getcwd()

# Read attributes for sample powder cans
flat_plate = pd.read_csv(os.path.join(this_dir, 'dict', 'flatPlate.csv'))
cylindrical = pd.read_csv(os.path.join(this_dir, 'dict', 'cylindrical.csv'))
annular = pd.read_csv(os.path.join(this_dir, 'dict', 'annular.csv'))

# Define remove_parentheses() function
def remove_parentheses(i):
    # Replace text within parentheses with empty string
    # \( matches an opening parenthesis
    # \) matches a closing parenthesis
    # [^)]* matches any character except a closing parenthesis
    return re.sub(r'\([^)]*\)', '', i) if i else None

# Define validate_formula() function
def validate_formula(format, formula):
    # Flag an error if formula (string) does not match format (pattern)
    if not re.fullmatch(format, formula):
        raise ValueError(f'The string does not match the accepted format.')
    # Return `True` if formula (string) matches format (pattern), `False` otherwise
    return True

# Define sum_total_n() function
def sum_total_n(i):
    # Define regular expression to find all digits in string
    digits = re.findall(r'\d', i)

    # Convert digit to integer and sum all digits
    total_sum = sum(int(j) for j in digits)
    
    # Return total sum
    return total_sum

# Define read_cif() function
def read_cif(filepath):
    # Read Crystallographic Information File
    cif = CifFile.ReadCif(filepath)

    # Define keys of interest
    keys = ['_chemical_formula_sum', '_cell_length_a', '_cell_length_b', '_cell_length_c', '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma', '_cell_volume', '_cell_formula_units_Z']

    # Initialize empty dictionary
    cif_dict = {}

    # Iterate over each block in CIF
    for database_code_PCD, block in cif.items():
        # Print Pearson's Crystal Data browser code
        print(f'Block: {database_code_PCD}')
    
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
                # Store `None` in dictionary
                block_data[i] = None
                print(f'{i}: Not found in {database_code_PCD}.')

        # Add block_data dictionary to information dictionary
        cif_dict[database_code_PCD] = block_data

    # Check for index troubleshooting
    # for database_code_PCD, block_data in cif_dict.items():
        # Print Pearson's Crystal Data browser code
        # print(f'Block: {database_code_PCD}')

    # Print value for each key of interest
    for key, value in block_data.items():
        print(f'{key}: {value}')
        print()

    # Extract variables of interest
    a = remove_parentheses(cif_dict[database_code_PCD]['_cell_length_a'])
    b = remove_parentheses(cif_dict[database_code_PCD]['_cell_length_b'])
    c = remove_parentheses(cif_dict[database_code_PCD]['_cell_length_c'])
    alpha = cif_dict[database_code_PCD]['_cell_angle_alpha']
    beta = cif_dict[database_code_PCD]['_cell_angle_beta']
    gamma = cif_dict[database_code_PCD]['_cell_angle_gamma']
    Z_param = cif_dict[database_code_PCD]['_cell_formula_units_Z']
    formula = cif_dict[database_code_PCD]['_chemical_formula_sum']
    cell_volume = cif_dict[database_code_PCD]['_cell_volume']

    # Find sample number density
    # Number density of sample in number of atoms or formula units per cubic Angstrom
    sample_n_density = eval(Z_param)/eval(cell_volume)

    # Split string by spaces
    formula_split = formula.split() if formula else []

    # Initialize empty dictionary
    elements = {}

    # Define regular expression pattern for dictionary
    # Group 1: [A-Za-z] matches any uppercase (A-Z) or lowercase (a-z) letter
    # Group 2: (\d+), where \d matches any digit (0-9) 
    pattern = re.compile(r'([A-Za-z]+)(\d+)')

    # Iterate through each element; return both index (i) and string for each element in chemical formula
    for i, string in enumerate(formula_split):
        # Match string against regular expression pattern
        match = pattern.match(string)
    
        # Check object is not `None`
        if match:
            # Extract element symbol/Group 1 match
            element_symbol = match.group(1)

            # Extract subscript/Group 2 match, if `None` (no subscript) then `1`
            # `None` is unlikely when retrieving from '_chemical_formula_sum' in CIF
            subscript = match.group(2) if match.group(2) else '1'

            # Find how many numbers per formula units
            n_fu = eval(Z_param)*eval(subscript)

            # Create new entry in dictionary for each element
            elements[f'element_{i+1}'] = {'symbol': element_symbol, 'subscript': subscript, 'n_fu': n_fu}

    # Print dictionary for troubleshooting
    # for key, value in elements.items():
        # print(f'{key}: {value}')

    # Initialize empty array
    array = []

    # Iterate over each item in dictionary
    for key, value in elements.items():
        # Retrieve element symbol and numbers per formula units 
        symbol = value['symbol']
        n_fu = value['n_fu']

        # Concatenate and append to array
        concatenate = f'{symbol}{n_fu}'
        array.append(concatenate)

    # Print array for troubleshooting
    # print(array)

    # Concatenate array to string
    mantid_formula = '-'.join([str(i) for i in array])

    # Print string for troubleshooting
    # print(mantid_formula)

    # Initialize variable
    total_n = 0

    # Iterate over dictionary items
    for key, value in elements.items():
        # Extract numbers per formula units value and add to total sum
        n_fu = value['n_fu']
        total_n += n_fu

    # Print total n
    # print('Total n:', total_n)
    
    # Return mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma
    return mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma

# Define xs_calculator function; cross-section calculator
def xs_calculator(x, neutron_energy, pack_fraction, can = 'flat', can_material = 'aluminum', Z_param = None, a = None, b = None, c = None, alpha = None, beta = None, gamma = None):
    # can = ['flat', 'cyl', 'annular']

    # Create empty data container/workspace
    ws = CreateSampleWorkspace()

    # Define pattern for validate_formula() function
    # Each element is followed by number of atoms for that element, specified without a hyphen; each element is separated from other elements using a hyphen
    # For isotopes, isotope must be enclosed by parenthesis; e.g., (Li7)
    mantid_format = r'^(\([A-Za-z][a-z]*\d*\)\d*|[A-Za-z][a-z]*\d*)(-[A-Za-z][a-z]*\d*|\([A-Za-z][a-z]*\d*\)\d*)*$'

    # Define valid file extension
    file_extension = '.cif'

    # Check if x is a filepath for a CIF
    try:
        if x.endswith(file_extension) and os.path.isfile(x):
            # Extract mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma with read_cif() function
            mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma = read_cif(x)
            
            # Define sample
            # Add material to data container/workspace
            SetSample(ws, Material = {'ChemicalFormula': mantid_formula,
                                      'SampleNumberDensity': sample_n_density})

    # If error, check if x is a string and validate with validate_formula() function
    except ValueError:
        if validate_formula(mantid_format, x):
            # Define mantid_formula if `True`
            mantid_formula = x

            # Define total_n if `True` with sum_total_n() function
            total_n = sum_total_n(mantid_formula)

            # Define sample
            # Add material to data container/workspace
            SetSample(ws, Material = {'ChemicalFormula': mantid_formula,
                                      'ZParameter': Z_param})
        else:
            print(f'{x} cannot be added to workspace. It does not match the accepted format nor is it a CIF.')

    # Check version for troubleshooting
    # mantid_version = mtd.__version__
    # print(mantid_version)

    # Check data container/workspace for troubleshooting
    # if ws.id() == 'Workspace2D':
    #    print(ws.name() + ' is an ' + ws.id())

    # Print formula for troubleshooting
    # print(mantid_formula)
    
    # Check for troubleshooting
    sample = ws.sample()

    # Retrieve absorption cross-section
    absorb_xs = sample.getMaterial().absorbXSection() * total_n
    # print(f'Absorption cross-section: {absorb_xs}')

    # Retrieve total scattering cross-section
    # Define scattering cross section per formula unit in bn/fu
    scatter_xs = sample.getMaterial().totalScatterXSection() * total_n
    # print(f'Total scattering cross-section: {scatter_xs}')

    # Retrieve coherent scattering length
    scatter_length = sample.getMaterial().cohScatterLength()
    # print(f'Coherent scattering length: {scatter_length}')

    # Retrieve relative molecular mass
    # Define molecular mass in g/mol/fu
    molecular_mass = sample.getMaterial().relativeMolecularMass()
    # print(f'Relative molecular mass: {molecular_mass}')

    # Convert string values to floats
    a = float(a)
    b = float(b)
    c = float(c)
    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)
    
    # Find unit cell volume in A^3
    unit_cell_volume = a * b * c * np.sqrt(1-np.cos(alpha * np.pi/180)**2 - np.cos(beta * np.pi/180)**2 - np.cos(gamma * np.pi/180)**2 + 2*np.cos(alpha * np.pi/180) * np.cos(beta * np.pi/180) * np.cos(gamma * np.pi/180))
    
    # Find theoretical density in g/cc
    theory_density = molecular_mass/unit_cell_volume/0.6022
    
    # Find absorption cross section per formula unit at given incident neutron energy in bn/fu
    sqrt_neutron_energy = np.sqrt(25/float(neutron_energy))
    absorb_xs_sqrt = absorb_xs*sqrt_neutron_energy
    
    # Initialize empty dictionaries
    sample_mass = {}
    sample_volume_dict = {}
    sample_moles = {}

    # Check if `flat` is in `can` parameter in xs_calculator() function
    if 'flat' in can:
        # Iterate over each row in DataFrame
        for index, row in flat_plate.iterrows():
            # Extract `drawing_number` and `sample_volume_value` for each row
            drawing_number = row['drawing_number']
            sample_volume_value = row['sample_volume']

            # Find sample mass
            sample_mass_flat = sample_volume_value*pack_fraction
        
            # Find sample volume in cc
            sample_volume_flat = sample_mass_flat/theory_density/pack_fraction

            # Find number of moles of formula unit in sample
            sample_moles_flat = sample_mass_flat/molecular_mass
            
            # Populate dictionaries with `drawing_number` as key
            sample_mass[drawing_number] = sample_mass_flat
            sample_volume_dict[drawing_number] = sample_volume_flat
            sample_moles[drawing_number] = sample_moles_flat

    # Check if `cyl` is in `can` parameter in xs_calculator() function
    # if 'cyl' in can:
        # Iterate over each row in DataFrame
        # for row in cylindrical.iterrows(): 
            # Find sample mass
            # sample_mass_cylind = 0

    # Check if `annular` is in `can` parameter in xs_calculator() function
    # if 'annular' in can:
        # Iterate over each row in DataFrame
        # for row in annular.iterrows():
            # Find sample mass
            # sample_mass_annular = sample_volume*pack_fraction
    
    # Find penetration depth due to scattering in cm
    scatter_depth = unit_cell_volume/(scatter_xs * pack_fraction)

    # if 'flat in can: 
    # Find thickness of a 10 percent scatterer in cm
    scatter_thick = np.log(0.9)*scatter_depth
    
    # Find penetration depth due to absorption in cm
    absorb_depth = unit_cell_volume/(absorb_xs_sqrt * pack_fraction)
    
    # Find total penetration depth due to scattering and absorption in cm
    total_depth = unit_cell_volume/((scatter_xs + absorb_xs_sqrt) * pack_fraction)

    # Initialize empty dictionaries
    sample_thick = {}
    percent_scatter = {}
    percent_absorb = {}

    # Check if `flat` is in `can` parameter in xs_calculator() function
    if 'flat' in can:
        # Iterate over each row in DataFrame
        for index, row in flat_plate.iterrows():
            # Extract `drawing_number` and `sample_area` for each row
            drawing_number = row['drawing_number']
            sample_area = row['sample_area']

            # Find thickness of sample spread homogenously over sample can in cm
            sample_thick_flat = sample_mass[drawing_number]/(pack_fraction*theory_density*sample_area) 
    
            # Find percent of incident beam that is scattered (assume no absorption)
            percent_scatter_flat = 100 * (1-(np.exp(-(sample_thick_flat/scatter_depth))))
    
            # Find percent of incident beam that is absorbed (assume no scattering)
            percent_absorb_flat =  100 * (1-(np.exp(-(sample_thick_flat/absorb_depth))))
            
            # Populate dictionaries with `drawing_number` as key
            sample_thick[drawing_number] = sample_thick_flat
            percent_scatter[drawing_number] = percent_scatter_flat
            percent_absorb[drawing_number] = percent_absorb_flat

    # Check if `cyl` is in `can` parameter in xs_calculator() function
    # if 'cyl' in can:
        # Find thickness of sample spread homogenously over sample can in cm
        # sample_thick_cyl = sample_mass_cyl/(pack_fraction*theory_density*can_volume) 
    
        # Find percent of incident beam that is scattered (assume no absorption)
        # percent_scatter_cyl = 100 * (1-(np.exp(-(sample_thick_cyl/scatter_depth))))
    
        # Find percent of incident beam that is absorbed (assume no scattering)
        # percent_absorb_cyl =  100 * (1-(np.exp(-(sample_thick_cyl/absorb_depth))))
    
    # Check if `annular` is in `can` parameter in xs_calculator() function
    # if 'annular' in can:
        # Find thickness of sample spread homogenously over sample can in cm
        # Find percent of incident beam that is scattered (assume no absorption)
        # Find percent of incident beam that is absorbed (assume no scattering)
 
    # Initialize empty dictionary to hold flat data
    flat_data = {
        'drawing_number': [],
        'sample_mass': [],
        'sample_volume': [],
        'sample_moles': [],
        'sample_thick': [],
        'percent_scatter': [],
        'percent_absorb': []
    }
    
    for id in drawing_number:
        flat_data['drawing_number'].append(id)
        flat_data['sample_mass'].append(sample_mass[id])
        flat_data['sample_volume'].append(sample_volume_dict[id])
        flat_data['sample_moles'].append(sample_moles[id])
        flat_data['sample_thick'].append(sample_thick[id])
        flat_data['percent_scatter'].append(percent_scatter[id])
        flat_data['percent_absorb'].append(percent_absorb[id])
 
    # Create a DataFrame
    df = pd.DataFrame(flat_data)
    print(df)

    # Print sample information
    print(f'''
    =============================
        Sample Information
    =============================
    Formula:          {mantid_formula}
    Z:                {Z_param:<10}
    Lattice constants:
    Packing fraction: {pack_fraction:<10}
    Incident neutron energy:   {neutron_energy:<10}
    ============================
    ''')

    # Print table headers
    # print(f'{headers_1[0]:<15} {headers_1[1]:<22} {headers_1[2]:<15} {headers_1[3]:<15} {headers_1[4]:<15} {headers_1[5]:<15} {headers_1[6]:<15}')

    # Print separator line
    # print('-' * 120)

# return scatter_depth, percent_scatter, absorb_depth, percent_absorb

cif_filepath = os.path.join(this_dir, 'cif', 'Ba2Co1La2O12Te2.cif')
xs_calculator(x = cif_filepath, neutron_energy = 5, pack_fraction = 0.45)
# xs_calculator(x = 'Ba2Co1La2O12Te2', neutron_energy = 5, pack_fraction = 0.45, Z_param = 3, a = 5.69, b = 5.69, c = 27.58, alpha = 90, beta = 90, gamma = 120)