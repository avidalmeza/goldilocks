# goldilocks.py

import os
import re
import numpy as np
import pandas as pd
import CifFile
import mantid as mtd
from mantid.simpleapi import CreateSampleWorkspace, SetSample
from mantid.kernel import Material
from scipy.integrate import quad

# Obtain current working directory filepath
this_dir = os.getcwd()

# Read attributes for sample powder cans
flat_plate = pd.read_csv(os.path.join(this_dir, 'src', 'dict', 'flatPlate.csv'))
cylindrical = pd.read_csv(os.path.join(this_dir, 'src', 'dict', 'cylindrical.csv'))
annulus = pd.read_csv(os.path.join(this_dir, 'src', 'dict', 'annulus.csv'))
cans_desc = pd.read_csv(os.path.join(this_dir, 'src', 'dict', 'canDescription.csv'))

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
    # Group 1: [A-Za-z]+ matches any uppercase (A-Z) or lowercase (a-z) letter or sequence of letters (i.e., element)
    # Group 2: (\([\w\d]+\)) matches any sequence denoted within parentheses (i.e., isotope like (Li7))
    # Group 3: (\d+(\.\d+)?) matches any digit (0-9), possibly with a decimal
    pattern = re.compile(r'([A-Za-z]+|\([\w\d]+\))(\d+(\.\d+)?)')

    # Iterate through each element; return both index (i) and string for each element in chemical formula
    for i, string in enumerate(formula_split):
        # Match string against regular expression pattern
        match = pattern.match(string)
    
        # Check object is not `None`
        if match:
            # Extract element/Group 1 match
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

    # Concatenate array to string and replace spaces with '-'
    mantid_formula = '-'.join([str(i) for i in array])
    # mantid_formula = mantid_formula.replace(' ', '-')

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
    return mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma, Z_param

# Define xs_calculator function; cross-section calculator
def xs_calculator(x, neutron_energy, pack_fraction, can = ['flat', 'cyl', 'annulus'], can_material = ['aluminum', 'vanadium'], Z_param = None, a = None, b = None, c = None, alpha = None, beta = None, gamma = None):
    # Create empty data container/workspace
    ws = CreateSampleWorkspace()

    # Define pattern for validate_formula() function
    # Each element is followed by number of atoms for that element, specified without a hyphen; each element is separated from other elements using a hyphen
    # For isotopes, isotope must be enclosed by parenthesis; e.g., (Li7)
    mantid_format = r'^(\([A-Za-z][a-z]*\d*\)\d*|[A-Za-z][a-z]*\d*)(-[A-Za-z][a-z]*\d*|\([A-Za-z][a-z]*\d*\)\d*)*$'

    # Define valid file extension
    file_extension = '.cif'

    # Initialize to avoid UnboundLocalError
    total_n = 0

    # Check if x is a filepath for a CIF
    try:
        if x.endswith(file_extension) and os.path.isfile(x):
            # Extract mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma with read_cif() function
            mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma, Z_param = read_cif(x)

            # Define sample in workspace
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

            # Convert string values to floats
            a = float(a)
            b = float(b)
            c = float(c)
            alpha = float(alpha)
            beta = float(beta)
            gamma = float(gamma)
    
            # Find unit cell volume in A^3
            unit_cell_volume = a * b * c * np.sqrt(1-np.cos(alpha * np.pi/180)**2 - np.cos(beta * np.pi/180)**2 - np.cos(gamma * np.pi/180)**2 + 2*np.cos(alpha * np.pi/180) * np.cos(beta * np.pi/180) * np.cos(gamma * np.pi/180))
    
            # Define sample
            # Add material to data container/workspace
            SetSample(ws, Material = {'ChemicalFormula': mantid_formula,
                                      'UnitCellVolume': float(unit_cell_volume),
                                      'ZParameter': float(Z_param)})
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
    
    # Obtain sample object from workspace
    sample = ws.sample()

    # Retrieve absorption cross-section
    # Define absorption cross section per formula unit in bn/fu
    absorb_xs = (float(sample.getMaterial().absorbXSection()*total_n))*np.sqrt(25/neutron_energy)
    # print(f'Absorption cross-section: {absorb_xs}')

    # Retrieve total scattering cross-section
    # Define scattering cross section per formula unit in bn/fu
    scatter_xs = float(sample.getMaterial().totalScatterXSection()*total_n)
    # print(f'Total scattering cross-section: {scatter_xs}')

    # Retrieve coherent scattering length
    scatter_length = float(sample.getMaterial().cohScatterLength())
    # print(f'Coherent scattering length: {scatter_length}')

    # Retrieve relative molecular mass
    # Define molecular mass in g/mol/fu
    molecular_mass = float(sample.getMaterial().relativeMolecularMass())
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

    # if 'alumnium' in can_material:
    #    element = 'Al'

    #if 'vanadium' in can_material:  
    #    element = 'V'
    
    # Initialize empty dictionaries
    sample_mass = {}
    can_volume_dict = {}
    sample_moles = {}

    # Check if `flat` is in `can` parameter in xs_calculator() function
    if 'flat' in can:
        # Iterate over each row in DataFrame
        for index, row in flat_plate.iterrows():
            # Extract `drawing_number` and `can_volume` for each row
            drawing_number = row['drawing_number']
            can_volume_flat = row['can_volume']/1000 # in cm^3

            # Find sample mass in grams
            sample_mass_flat = can_volume_flat * theory_density * pack_fraction # in grams
        
            # Find number of moles of formula unit in sample
            sample_moles_flat = (sample_mass_flat/molecular_mass)*1000
            
            # Populate dictionaries with `drawing_number` as key
            sample_mass[drawing_number] = sample_mass_flat
            can_volume_dict[drawing_number] = can_volume_flat
            sample_moles[drawing_number] = sample_moles_flat

    # Check if `cyl` is in `can` parameter in xs_calculator() function
    if 'cyl' in can:
    # Iterate over each row in DataFrame
        for index, row in cylindrical.iterrows():
            # Extract `drawing_number` and `sample_volume` for each row
            drawing_number = row['drawing_number']
            can_volume_cyl = row['can_volume']

            # Find sample mass in grams
            sample_mass_cyl = can_volume_cyl * theory_density * pack_fraction # in grams

            # Find number of moles of formula unit in sample
            sample_moles_cyl = (sample_mass_cyl/molecular_mass)*1000

            # Populate dictionaries with `drawing_number` as key
            sample_mass[drawing_number] = sample_mass_cyl
            sample_moles[drawing_number] = sample_moles_cyl
            can_volume_dict[drawing_number] = can_volume_cyl

    # Check if `annulus` is in `can` parameter in xs_calculator() function
    if 'annulus' in can:
        # Iterate over each row in DataFrame
        for index, row in annulus.iterrows():
            # Extract `drawing_number` and `sample_volume` for each row
            drawing_number = row['drawing_number']
            can_volume_ann = row['sample_volume'] # Note: update name in dict

            # Find sample mass in grams
            sample_mass_ann = can_volume_ann * theory_density * pack_fraction # in grams

            # Find number of moles of formula unit in sample
            sample_moles_ann = (sample_mass_cyl/molecular_mass)*1000

            # Populate dictionaries with `drawing_number` as key
            sample_mass[drawing_number] = sample_mass_ann
            sample_moles[drawing_number] = sample_moles_ann
            can_volume_dict[drawing_number] = can_volume_ann
    
    # Find penetration depth due to scattering in cm
    scatter_depth = unit_cell_volume/(scatter_xs * pack_fraction)

    # Find thickness of a 10 percent scatterer in cm
    scatter_thick = np.log(0.9)*scatter_depth
    
    # Find penetration depth due to absorption in cm
    absorb_depth = unit_cell_volume/(absorb_xs * pack_fraction)
    
    # Find total penetration depth due to scattering and absorption in cm
    total_depth = unit_cell_volume/((scatter_xs + absorb_xs) * pack_fraction)

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

            # Find thickness of sample spread homogenously over sample can in mm
            sample_thick_flat = row['sample_thick']/10

            # Extract sample mass from dictionary
            sample_mass_i = sample_mass.get(drawing_number)

            # Find percent of incident beam that is scattered (assume no absorption)
            percent_scatter_flat = 100 * (1-(np.exp(-(sample_thick_flat/scatter_depth))))
    
            # Find percent of incident beam that is absorbed (assume no scattering)
            percent_absorb_flat =  100 * (1-(np.exp(-(sample_thick_flat/absorb_depth))))
            
            # Populate dictionaries with `drawing_number` as key
            sample_thick[drawing_number] = sample_thick_flat
            percent_scatter[drawing_number] = percent_scatter_flat
            percent_absorb[drawing_number] = percent_absorb_flat

    # Check if `cyl` is in `can` parameter in xs_calculator() function
    if 'cyl' in can:
        # Iterate over each row in DataFrame
        for index, row in cylindrical.iterrows():
            # Extract `drawing_number`, `sample_area`, and `can_thick` for each row
            drawing_number = row['drawing_number']
            can_inner_radius = row['can_inner_radius']
            can_thick = row['can_thick']

            # Define integrand_cyl() function
            def integrand_cyl(x, paramwave):
                R1 = paramwave[0]
                zeta = paramwave[1]

                dval = np.sqrt(R1**2 - x**2)
                integrand = 1 - np.exp(-2 * dval / zeta)

                # Return integrand for cylindrical integration
                return integrand

            # Define integral_cyl() function
            def integral_cyl(xs, unit_cell_volume, pack_fraction, R1):
                zeta = unit_cell_volume/pack_fraction/xs # Note: packing efficiency? units?
                paramwave = [R1, zeta]
                int1, _ = quad(integrand_cyl, -R1, R1, args = (paramwave))
                result = 100 * int1/(2 * R1)
                return result
            
            # Find percent of incident beam that is scattered (assume no absorption)
            percent_scatter_cyl = integral_cyl(xs = scatter_depth, unit_cell_volume = unit_cell_volume, pack_fraction = pack_fraction, R1 = can_inner_radius)
            
            # Find percent of incident beam that is absorbed (assume no scattering)
            percent_absorb_cyl = integral_cyl(xs = absorb_depth, unit_cell_volume = unit_cell_volume, pack_fraction = pack_fraction, R1 = can_inner_radius)

            # Populate dictionaries with `drawing_number` as key
            percent_scatter[drawing_number] = percent_scatter_cyl
            percent_absorb[drawing_number] = percent_absorb_cyl

    # Check if `annulus` is in `can` parameter in xs_calculator() function
    # if 'annulus' in can:
        # Define integrand_ann() function
        # def integrand_ann(inx, pwave):
        #     R1 = pwave[0]
        #     zeta = pwave[1]
        #     R2 = pwave[2]
            
        #     if abs(inx) < R1:
        #         return 0
        #     elif abs(inx) > R2:
        #         return 0
            
        #     dval = np.sqrt(R2**2 - inx**2)

        #     integrand = 1 - np.exp(-2 * dval / zeta)
            
        #     return integrand

        # # Define integral_ann() function
        # def integral_ann(xs, unit_cell_volume, pack_fraction, R1, R2):
        #     zeta = unit_cell_volume/pack_fraction/xs
        #     paramwave = [R1, zeta, R2]
            
        #     int1, _ = quad(integrand_ann, -R2, R2, args = (paramwave))
        #     result = 100 * int1 / (2 * R2)
            
        #     return result
        # Find percent of incident beam that is scattered (assume no absorption)
        # Find percent of incident beam that is absorbed (assume no scattering)

    # Convert dictionaries to DataFrames
    # Create dictionary `df_dict` of DataFrames
    df_dict = {
        'sample_mass': pd.DataFrame(sample_mass.items(), columns = ['drawing_number', 'sample_mass']),
        'can_volume': pd.DataFrame(can_volume_dict.items(), columns = ['drawing_number', 'can_volume']),
        'sample_moles': pd.DataFrame(sample_moles.items(), columns = ['drawing_number', 'sample_moles']),
        'sample_thick': pd.DataFrame(sample_thick.items(), columns = ['drawing_number', 'sample_thick']),
        'percent_scatter': pd.DataFrame(percent_scatter.items(), columns = ['drawing_number', 'percent_scatter']),
        'percent_absorb': pd.DataFrame(percent_absorb.items(), columns = ['drawing_number', 'percent_absorb'])
    }

    # Extract DataFrame `sample_mass` from dictionary
    df = df_dict['sample_mass']

    # Merge with left join on `drawing_number`
    for key in list(df_dict.keys())[1:]:
        df = df.merge(df_dict[key], on = 'drawing_number', how = 'left')

    # Merge with left join on `drawing_number`
    df_concat = pd.merge(cans_desc, df, on = 'drawing_number', how = 'left')

    # Drop `drawing_number` column
    df_concat.drop('drawing_number', axis = 1, inplace = True)

    # Print sample information
    print(f'''
    =============================
        Sample Information
    =============================
    Formula:          {mantid_formula}
    Z:                {Z_param}
    Lattice constant (a): {a} Å
    Lattice constant (b): {b} Å
    Lattice constant (c): {c} Å
    Lattice constant (alpha): {alpha}°
    Lattice constant (beta): {beta}°
    Lattice constant (gamma): {gamma}°
    Packing fraction: {pack_fraction}
    Incident neutron energy:   {neutron_energy} meV
    ============================
    ''')

    # Print sample can independent information
    print(f'''
    =========================================
        Sample Can Independent Values
    =========================================
    Calculated unit cell volume (\Angs^3): {round(unit_cell_volume, 3)}
    Scattering cross section (bn/fu): {round(scatter_xs, 3)}
    Absorption cross section (bn/fu): {round(absorb_xs, 3)}
    Scattering penetration depth (cm): {round(scatter_depth, 3)}
    Absorption penetration depth (cm): {round(absorb_depth, 3)}
    Total peneteration depth (cm):  {round(total_depth, 3)}
    Calculated sample density (g/cc): {round(theory_density, 3)}
    ==========================================
    ''')

    print(df_concat)