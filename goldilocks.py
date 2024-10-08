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
from scipy.integrate import trapezoid
import datetime
import csv

# Obtain current working directory filepath
this_dir = os.getcwd()

# Read dictionaries as DataFrame
flat_plate = pd.read_csv(os.path.join(this_dir, "src", "dict", "flatPlate.csv"))
cylindrical = pd.read_csv(os.path.join(this_dir, "src", "dict", "cylindrical.csv"))
annulus = pd.read_csv(os.path.join(this_dir, "src", "dict", "annulus.csv"))
cans_desc = pd.read_csv(os.path.join(this_dir, "src", "dict", "canDescription.csv"))
material = pd.read_csv(os.path.join(this_dir, "src", "dict", "material.csv"))

# Define remove_parentheses() function
def remove_parentheses(i):
    # Replace text within parentheses (including parentheses) with empty string
    # \( matches an opening parenthesis
    # \) matches a closing parenthesis
    # [^)]* matches any character except a closing parenthesis
    return re.sub(r"\([^)]*\)", "", i) if i else None

# Define validate_formula() function
def validate_formula(format, formula):
    # Flag an error if formula (string) does not match format (pattern)
    if not re.fullmatch(format, formula):
        raise ValueError(f"The string does not match the accepted format.")
    # Return True if formula (string) matches format (pattern), False otherwise
    return True

# Define sum_total_n() function, used for molecular formula case
def sum_total_n(s):
    # Remove substrings within parentheses (i.e., isotopes) and extract remaining numeric values
    s = re.sub(r"\(.*?\)", "", s)
    numbers = re.findall(r"\d+\.?\d*", s)

    # Convert each number found to a float and sum
    total_sum = sum(float(num) for num in numbers)
    return total_sum

# Define uc_volume() function
def uc_volume(a, b, c, alpha, gamma, beta):
    # Find unit cell volume in A^3
    unit_cell_volume = a * b * c * np.sqrt(1 - np.cos(alpha * np.pi / 180)**2 - np.cos(beta * np.pi / 180)**2 - np.cos(gamma * np.pi / 180)**2 + 2*np.cos(alpha * np.pi / 180) * np.cos(beta * np.pi / 180) * np.cos(gamma * np.pi / 180))
    return unit_cell_volume

# Define material_properties() function
def material_properties(absorb_xs, scatter_xs, molecular_mass, unit_cell_volume, Z_param):
    # Find penetration depth due to scattering in cm
    material_scatter_depth = unit_cell_volume / scatter_xs # in cm
    
    # Find penetration depth due to absorption in cm
    material_absorb_depth = unit_cell_volume / absorb_xs # in cm
    
    # Find total penetration depth due to scattering and absorption in cm
    material_total_depth = unit_cell_volume / (scatter_xs + absorb_xs) # in cm
    
    # Find thickness of a 10 percent scatterer in cm
    material_scatter_thick = np.log(0.9) * material_scatter_depth # in cm
    
    # Find theoretical density in g/cc
    material_theory_density = (molecular_mass / unit_cell_volume / 0.6022) * Z_param # in g/cc
    
    # Return material_scatter_depth, material_absorb_depth, material_total_depth, material_scatter_thick, material_theory_density
    return material_scatter_depth, material_absorb_depth, material_total_depth, material_scatter_thick, material_theory_density

# Define integrand_cyl() function
def integrand_cyl(x, paramwave):
    # Extract inner can radius (R1) and zeta
    R1 = paramwave[0]
    zeta = paramwave[1]

    dval = np.sqrt(R1**2 - x**2)
    integrand = 1 - np.exp(-2 * dval / zeta)
    return integrand

# Define integral_cyl() function
def integral_cyl(xs, unit_cell_volume, pack_fraction, R1):
    zeta = unit_cell_volume/pack_fraction/xs # Note: packing efficiency? units?
    paramwave = [R1, zeta]
    int1, _ = quad(integrand_cyl, -R1, R1, args = (paramwave))
    result = 100 * int1 / (2 * R1)
    return result

# Define integrand_ann() function
def integrand_ann(x, paramwave):
    # Extract inner can radius (R1), zeta, and outer can radius (R2)
    R1 = paramwave[0]
    zeta = paramwave[1]
    R2 = paramwave[2]

    # If absolute value of x is less than inner can radius
    if np.abs(x) < R1:
        dval1 = np.sqrt(R1**2 - x**2)
    # If absolute value of x is greater than outer can radius
    else:
        dval1 = 0
    
    dval2 = np.sqrt(R2**2 - x**2)
    integrand = 1 - np.exp(-2 * (dval2 - dval1) / zeta)
    return integrand

# Define integral_ann() function
def integral_ann(xs, unit_cell_volume, pack_fraction, R1, R2):
    zeta = unit_cell_volume / pack_fraction / xs
    paramwave = [R1, zeta, R2]
    
    # Set number of points
    n_points = 1000
    # Set integration limits
    x_values = np.linspace(-R2, R2, n_points)

    # Evaluate integrand at each point
    integrand_values = np.array([integrand_ann(x, paramwave) for x in x_values])
    
    # Perform trapezoidal integration
    integration = 100 * trapezoid(integrand_values, x_values)
    result = integration / (2 * R2)
    return result            

# Define calculate_flatPlate() function
def calculate_flatPlate(can_total_thickness, can_volume, scatter_depth, absorb_depth, theory_density):
    # Find percent of incident beam that is scattered (assume no absorption)
    can_percent_scatter = 100 * (1-(np.exp(-(can_total_thickness / scatter_depth))))

    # Find percent of incident beam that is absorbed (assume no scattering)
    can_percent_absorb = 100 * (1-(np.exp(-(can_total_thickness / absorb_depth))))

    # Find sample mass in grams
    can_mass = can_volume * theory_density # in grams

    # Return can_percent_scatter, can_percent_absorb, can_mass
    return can_percent_scatter, can_percent_absorb, can_mass

# Define can_annulus() function
def can_annulus(scatter_xs, absorb_xs, theory_density, unit_cell_volume, height, inner_radius, outer_radius, volume):
    # Find percent of incident beam that is scattered (assume no absorption)
    can_percent_scatter = integral_ann(scatter_xs, unit_cell_volume, pack_fraction = 1.0, R1 = inner_radius, R2 = outer_radius)
    
    # Find percent of incident beam that is absorbed (assume no scattering)
    can_percent_absorb = integral_ann(absorb_xs, unit_cell_volume, pack_fraction = 1.0, R1 = inner_radius, R2 = outer_radius)
    
    # Find sample mass in grams
    volume = height * np.pi * (outer_radius**2 - inner_radius**2) # in cm^3
    can_mass = volume * theory_density # in grams

    # Return can_percent_scatter, can_percent_absorb, can_mass
    return can_percent_scatter, can_percent_absorb, can_mass

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
    #    print(f"{key}: {value}")
    #    print()

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
    sample_n_density = eval(Z_param)/eval(cell_volume)

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
            # None is unlikely when retrieving from _chemical_formula_sum in CIF
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

# Define xs_calculator function; cross-section calculator
def xs_calculator(x, pack_fraction = 1.0, neutron_energy = None, neutron_wavelength = None, can = ["flat", "cyl", "annulus"], density = None, Z_param = None, a = None, b = None, c = None, alpha = None, beta = None, gamma = None):
    # Create empty data container/workspace
    ws = CreateSampleWorkspace()

    # Define pattern for validate_formula() function
    # Each element is followed by number of atoms for that element, specified without a hyphen; each element is separated from other elements using a hyphen
    # For isotopes, isotope must be enclosed by parenthesis; e.g., (Li7)
    # [A-Za-z][a-z]* matches a capital letter (e.g., C, O, B, P, F) or a capital letter followed by lowercase letters (e.g., Cu, Si, Al)
    # \d*(\.\d+)? matches any digit (0-9), possibly with a decimal
    # \([A-Za-z][a-z]?\d*(\.\d+)?\)\d*(\.\d+)? matches an isotope
    mantid_format = r"^([A-Za-z][a-z]?\d*(\.\d+)?|\([A-Za-z][a-z]?\d*(\.\d+)?\)\d*(\.\d+)?)(-([A-Za-z][a-z]?\d*(\.\d+)?|\([A-Za-z][a-z]?\d*(\.\d+)?\)\d*(\.\d+)?))*$"

    # Define valid file extension
    file_extension = ".cif"

    # Initialize to avoid UnboundLocalError
    mantid_formula = "" # Initialize type str
    total_n = 0

    # Check if x is a filepath for a CIF
    try:
        if x.endswith(file_extension) and os.path.isfile(x):
            # Gather mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma from read_cif() function
            mantid_formula, sample_n_density, total_n, a, b, c, alpha, beta, gamma, Z_param = read_cif(x)

            # Define sample in workspace
            # Add material to data container/workspace
            SetSample(ws, Material = {"ChemicalFormula": mantid_formula,
                                      "SampleNumberDensity": sample_n_density})
        else:
            raise ValueError("File is not a CIF or does not exist.")


    # If error, check if x is a string and validate with validate_formula() function
    except (ValueError, AttributeError) as e:
        # If isinstance() and validate_formula() returns True 
        if isinstance(x, str) and validate_formula(mantid_format, x):
            # Define mantid_formula if True
            mantid_formula = str(x)

            # Define total_n if True with sum_total_n() function
            total_n = sum_total_n(mantid_formula)
            
            # Check density is not None
            if density is not None:
                # Set false Z_param
                Z_param = 1

                # Define sample
                # Add material to data container/workspace
                SetSample(ws, Material = {"ChemicalFormula": str(mantid_formula),
                                          "SampleMassDensity": float(density)}) # in g/cm^3
                
                # Obtain sample object from workspace
                sample = ws.sample()

                # Retrieve relative molecular mass
                # Define molecular mass in g/mol/fu
                molecular_mass = float(sample.getMaterial().relativeMolecularMass()) # in g/mol/fu

                # Set false lattice constants, assume cubic structure
                a = np.power(molecular_mass / density / 0.6022, 1 / 3)
                b = np.power(molecular_mass / density / 0.6022, 1 / 3)
                c = np.power(molecular_mass / density / 0.6022, 1 / 3)
                alpha = 90
                beta = 90
                gamma = 90

                # Find unit cell volume in A^3
                unit_cell_volume = uc_volume(a, b, c, alpha, beta, gamma)
            else:
                # Ensure parameters are provided, raise an error if any are None
                if any(param is None for param in [a, b, c, alpha, beta, gamma, Z_param]):
                    raise ValueError("Parameters a, b, c, alpha, beta, gamma, and Z_param required if x is a string.")

                # Convert string values to floats
                a = float(a)
                b = float(b)
                c = float(c)
                alpha = float(alpha)
                beta = float(beta)
                gamma = float(gamma)

                # Find unit cell volume in A^3
                unit_cell_volume = uc_volume(a, b, c, alpha, beta, gamma)

                # Define sample
                # Add material to data container/workspace
                SetSample(ws, Material = {"ChemicalFormula": str(mantid_formula),
                                          "UnitCellVolume": float(unit_cell_volume),
                                          "ZParameter": float(Z_param)})
        else:
            print(f"{x} cannot be added to workspace. It does not match the accepted format nor is it a CIF.")
            # Exit function if x is invalid
            return 

    # Check version for troubleshooting
    # mantid_version = mtd.__version__
    # print(mantid_version)

    # Check data container/workspace for troubleshooting
    # if ws.id() == "Workspace2D":
    #    print(ws.name() + " is an " + ws.id())

    # Print formula for troubleshooting
    # print(mantid_formula)

    # Initialize to avoid UnboundLocalError
    neutron_txt = ""

    # If neutron_wavelength is given, calculate neutron_energy 
    if neutron_energy is None and neutron_wavelength is not None:
        neutron_energy = 81.81 / neutron_wavelength / neutron_wavelength # in meV

    # If neutron_energy is given, calculate neutron_wavelength 
    if neutron_energy is not None and neutron_wavelength is None:
        neutron_wavelength = np.sqrt(81.81 / neutron_energy) # in Å

    # If neutron_energy and neutron_wavelength are not given
    if neutron_energy is None and neutron_wavelength is None:
        neutron_energy = 14.7 # in meV
        neutron_wavelength = 2.359 # in Å

        neutron_txt = f"""
        A value for neutron energy or wavelength is not given. Neutron energy set to 14.7 meV and neutron wavelength set to 2.359 Angstroms.
        """
    
    """
    Find scattering and absorption for sample
    """
    # Obtain sample object from workspace
    sample = ws.sample()

    # Retrieve absorption cross-section
    # Define absorption cross section per formula unit in bn/fu
    absorb_xs = (float(sample.getMaterial().absorbXSection() * total_n)) * np.sqrt(25 / neutron_energy) # in bn/fu
    # print(f"Absorption cross-section: {absorb_xs}")

    # Retrieve total scattering cross-section
    # Define scattering cross section per formula unit in bn/fu
    scatter_xs = float(sample.getMaterial().totalScatterXSection() * total_n) # in bn/fu
    # print(f"Total scattering cross-section: {scatter_xs}")

    # Retrieve incoherent scattering cross-section
    # Define incoherent scattering cross section per formula unit in bn/fu
    incoh_scatter_xs = float(sample.getMaterial().incohScatterXSection() * total_n) # in bn/fu

    # Retrieve coherent scattering cross-section
    # Define coherent scattering cross section per formula unit in bn/fu
    coh_scatter_xs = float(sample.getMaterial().cohScatterXSection() * total_n) # in bn/fu

    # Retrieve coherent scattering length
    scatter_length = float(sample.getMaterial().cohScatterLength()) # in fm
    # print(f"Coherent scattering length: {scatter_length}")

    # Retrieve relative molecular mass
    # Define molecular mass in g/mol/fu
    molecular_mass = float(sample.getMaterial().relativeMolecularMass()) # in g/mol/fu
    # print(f"Relative molecular mass: {molecular_mass}")

    # Find equivelent mass of vanadium per mass of sample
    molecular_mass_vanadium = 50.941 # in g/mol/atom of vanadium
    incoh_xs_vanadium = 5.08 # in barns/atom 
    van_equiv_incoherent = molecular_mass_vanadium / incoh_xs_vanadium * incoh_scatter_xs / molecular_mass

    # Convert string values to floats
    a = float(a)
    b = float(b)
    c = float(c)
    alpha = float(alpha)
    beta = float(beta)
    gamma = float(gamma)
    
    # Find unit cell volume in A^3
    unit_cell_volume = uc_volume(a, b, c, alpha, beta, gamma)
    
    # Find theoretical density in g/cc
    theory_density = molecular_mass / unit_cell_volume / 0.6022 # in g/cc

    # Determine neutron wavelength based upon neutron energy
    neutron_wavelength = np.sqrt(81.81 / neutron_energy)

    # Initialize empty dictionaries
    sample_mass = {}
    can_volume_dict = {}
    sample_moles = {}

    # Check if flat is in can parameter in xs_calculator() function
    if "flat" in can:
        # Iterate over each row in DataFrame
        for index, row in flat_plate.iterrows():
            # Extract id and can_volume_mm3 for each row
            id = row["id"]
            can_volume_flat = row["can_volume_mm3"] / 1000 # in cm^3

            # Find sample mass in grams
            sample_mass_flat = can_volume_flat * theory_density * pack_fraction # in grams
        
            # Find number of moles of formula unit in sample
            sample_moles_flat = sample_mass_flat / molecular_mass

            # Populate dictionaries with id as key
            sample_mass[id] = round(sample_mass_flat, 3)
            can_volume_dict[id] = round(can_volume_flat, 3)
            sample_moles[id] = sample_moles_flat

    # Check if cyl is in can parameter in xs_calculator() function
    if "cyl" in can:
    # Iterate over each row in DataFrame
        for index, row in cylindrical.iterrows():
            # Extract id and can_volume_mm3 for each row
            id = row["id"]
            can_volume_cyl = row["can_volume_mm3"] / 1000 # in cm^3

            # Find sample mass in grams
            sample_mass_cyl = can_volume_cyl * theory_density * pack_fraction # in grams

            # Find number of moles of formula unit in sample
            sample_moles_cyl = sample_mass_cyl / molecular_mass

            # Populate dictionaries with id as key
            sample_mass[id] = round(sample_mass_cyl, 3)
            sample_moles[id] = sample_moles_cyl
            can_volume_dict[id] = round(can_volume_cyl, 3)

    # Check if annulus is in can parameter in xs_calculator() function
    if "annulus" in can:
        # Iterate over each row in DataFrame
        for index, row in annulus.iterrows():
            # Extract id and can_volume_mm3 for each row
            id = row["id"]
            can_volume_ann = row["sample_volume_mm3"] / 1000 # in cm^3

            # Find sample mass in grams
            sample_mass_ann = can_volume_ann * theory_density * pack_fraction # in grams

            # Find number of moles of formula unit in sample
            sample_moles_ann = sample_mass_cyl / molecular_mass

            # Populate dictionaries with id as key
            sample_mass[id] = round(sample_mass_ann, 3)
            sample_moles[id] = sample_moles_ann
            can_volume_dict[id] = round(can_volume_ann, 3)
    
    # Find penetration depth due to scattering in cm
    scatter_depth = unit_cell_volume / (scatter_xs * pack_fraction)

    # Find thickness of a 10 percent scatterer in cm
    scatter_thick = np.log(0.9) * scatter_depth
    
    # Find penetration depth due to absorption in cm
    absorb_depth = unit_cell_volume / (absorb_xs * pack_fraction)
    
    # Find total penetration depth due to scattering and absorption in cm
    total_depth = unit_cell_volume / ((scatter_xs + absorb_xs) * pack_fraction)

    # Initialize empty dictionaries
    percent_scatter = {}
    percent_absorb = {}

    # Check if flat is in can parameter in xs_calculator() function
    if "flat" in can:
        # Iterate over each row in DataFrame
        for index, row in flat_plate.iterrows():
            # Extract id for each row
            id = row["id"]

            # Find thickness of sample spread homogenously over sample can
            sample_thick_flat = row["sample_thick_mm"] / 10 # in cm

            # Find percent of incident beam that is scattered (assume no absorption)
            percent_scatter_flat = 100 * (1-(np.exp(-(sample_thick_flat / scatter_depth))))
    
            # Find percent of incident beam that is absorbed (assume no scattering)
            percent_absorb_flat =  100 * (1-(np.exp(-(sample_thick_flat / absorb_depth))))
            
            # Populate dictionaries with id as key
            percent_scatter[id] = round(percent_scatter_flat, 2)
            percent_absorb[id] = round(percent_absorb_flat, 2)

    # Check if cyl is in can parameter in xs_calculator() function
    if "cyl" in can:
        # Iterate over each row in DataFrame
        for index, row in cylindrical.iterrows():
            # Extract id and can_inner_radius_mm for each row
            id = row["id"]
            can_inner_radius = row["can_inner_radius_mm"] / 10 # in cm

            # Find percent of incident beam that is scattered (assume no absorption)
            percent_scatter_cyl = integral_cyl(scatter_xs, unit_cell_volume, pack_fraction, R1 = can_inner_radius)
            
            # Find percent of incident beam that is absorbed (assume no scattering)
            percent_absorb_cyl = integral_cyl(absorb_xs, unit_cell_volume, pack_fraction, R1 = can_inner_radius)

            # Populate dictionaries with id as key
            percent_scatter[id] = round(percent_scatter_cyl, 2)
            percent_absorb[id] = round(percent_absorb_cyl, 2)

    # Check if annulus is in can parameter in xs_calculator() function
    if "annulus" in can:
        # Iterate over each row in DataFrame
        for index, row in annulus.iterrows():
            # Extract id, sample_inner_radius_mm, and sample_outer_radius_mm for each row
            id = row["id"]
            sample_inner_radius_ann = row["sample_inner_radius_mm"] / 10 # in cm
            sample_outer_radius_ann = row["sample_outer_radius_mm"] / 10 # in cm

            # Find percent of incident beam that is scattered (assume no absorption)
            percent_scatter_ann = integral_ann(scatter_xs, unit_cell_volume, pack_fraction, R1 = sample_inner_radius_ann, R2 = sample_outer_radius_ann)

            # Find percent of incident beam that is absorbed (assume no scattering)
            percent_absorb_ann = integral_ann(absorb_xs, unit_cell_volume, pack_fraction, R1 = sample_inner_radius_ann, R2 = sample_outer_radius_ann)

            # Populate dictionaries with id as key
            percent_scatter[id] = round(percent_scatter_ann, 2)
            percent_absorb[id] = round(percent_absorb_ann, 2)

    """
    Find scattering and absorption for sample can
    """
    # Initialize empty dictionaries
    can_percent_scatter = {}
    can_percent_absorb = {}
    can_mass = {}

    # Check if flat is in can parameter in xs_calculator() function
    if "flat" in can:
        # Iterate over each row in DataFrame
        for index, row in flat_plate.iterrows():
            # Extract id, material, can_total_thick_mm, and can_volume_mm3 for each row
            id = row["id"]
            can_material = row["material"]
            can_total_thickness = row["can_total_thick_mm"] / 10  # in cm
            can_volume = row["can_material_volume_cm3"]

            material_row = material[material["material"] == can_material].iloc[0]
            
            # Set absorption cross-section in bn/fu
            material_absorb_xs = material_row["absorb_xs"] * np.sqrt(25 / neutron_energy) * material_row["Z_param"] # in bn/fu
            # Set total scattering cross-section in bn/fu
            material_scatter_xs = material_row["scatter_xs"] * material_row["Z_param"] # in bn/fu
            # Set relative molecular mass
            material_molecular_mass = material_row["molecular_mass"]
            # Set unit cell volume
            material_unit_cell_volume = material_row["unit_cell_volume"]

            # Gather material_scatter_depth, material_absorb_depth, material_theory_density from material_properties() function
            material_scatter_depth, material_absorb_depth, _, _, material_theory_density = material_properties(material_absorb_xs, material_scatter_xs, material_molecular_mass, material_unit_cell_volume, material_row["Z_param"])
            
            # Find percent of incident beam that is scattered (assume no absorption)
            # Find percent of incident beam that is absorbed (assume no scattering)
            # Find can mass in grams
            can_percent_scatter_flat, can_percent_absorb_flat, can_mass_flat = calculate_flatPlate(can_total_thickness = can_total_thickness, can_volume = can_volume, scatter_depth = material_scatter_depth, absorb_depth = material_absorb_depth, theory_density = material_theory_density)

            # Populate dictionaries with id as key
            can_percent_scatter[id] = round(can_percent_scatter_flat, 2)
            can_percent_absorb[id] = round(can_percent_absorb_flat, 2)
            can_mass[id] = round(can_mass_flat, 3)

    # Check if cyl is in can parameter in xs_calculator() function
    if "cyl" in can:
        # Iterate over each row in DataFrame
        for index, row in cylindrical.iterrows():
            # Extract id, material, can_inner_radius_mm, can_outer_radius_mm, and sample_height_mm for each row
            id = row["id"]
            can_material = row["material"]
            can_inner_radius = row["can_inner_radius_mm"] / 10 # in cm
            can_outer_radius = row["can_outer_radius_mm"] / 10 # in cm
            sample_height = row["sample_height_mm"] / 10 # in cm
            can_material_volume = row["can_material_volume_cm3"]

            material_row = material[material["material"] == can_material].iloc[0]
            
            # Set absorption cross-section in bn/fu
            material_absorb_xs = material_row["absorb_xs"] * np.sqrt(25 / neutron_energy) * material_row["Z_param"] # in bn/fu
            # Set total scattering cross-section in bn/fu
            material_scatter_xs = material_row["scatter_xs"] * material_row["Z_param"] # in bn/fu
            # Set relative molecular mass
            material_molecular_mass = material_row["molecular_mass"]
            # Set unit cell volume
            material_unit_cell_volume = material_row["unit_cell_volume"]
            
            # Gather material_scatter_depth, material_absorb_depth, material_theory_density from material_properties() function
            material_scatter_depth, material_absorb_depth, _, _, material_theory_density = material_properties(material_absorb_xs, material_scatter_xs, material_molecular_mass, material_unit_cell_volume, material_row["Z_param"])
            
            # Find percent of incident beam that is scattered (assume no absorption)
            # Find percent of incident beam that is absorbed (assume no scattering)
            # Find can mass in grams
            can_percent_scatter_cyl, can_percent_absorb_cyl, can_mass_cyl = can_annulus(scatter_xs = material_scatter_xs, absorb_xs = material_absorb_xs, theory_density = material_theory_density, unit_cell_volume = material_unit_cell_volume, height = sample_height, inner_radius = can_inner_radius, outer_radius = can_outer_radius, volume = can_material_volume)

            # Populate dictionaries with id as key
            can_percent_scatter[id] = round(can_percent_scatter_cyl, 2)
            can_percent_absorb[id] = round(can_percent_absorb_cyl, 2)
            can_mass[id] = round(can_mass_cyl, 3)

    # Check if annulus is in can parameter in xs_calculator() function
    if "annulus" in can:
        # Iterate over each row in DataFrame
        for index, row in annulus.iterrows():
            # Extract id, material, can_R3_cm, new_can_R4_cm, can_height_cm, and new_can_material_total_volume_cm3 for each row
            id = row["id"]
            can_material = row["material"]
            can_R3_cm = row["can_R3_cm"]
            new_can_R4_cm = row["new_can_R4_cm"]
            can_height_cm = row["can_height_cm"]
            can_material_volume_cm3 = row["new_can_material_total_volume_cm3"]

            material_row = material[material["material"] == can_material].iloc[0]
            
            # Set absorption cross-section in bn/fu
            material_absorb_xs = material_row["absorb_xs"] * np.sqrt(25 / neutron_energy) * material_row["Z_param"] # in bn/fu
            # Set total scattering cross-section in bn/fu
            material_scatter_xs = material_row["scatter_xs"] * material_row["Z_param"] # in bn/fu
            # Set relative molecular mass
            material_molecular_mass = material_row["molecular_mass"]
            # Set unit cell volume
            material_unit_cell_volume = material_row["unit_cell_volume"]

            # Gather material_scatter_depth, material_absorb_depth, material_theory_density from material_properties() function
            material_scatter_depth, material_absorb_depth, _, _, material_theory_density = material_properties(material_absorb_xs, material_scatter_xs, material_molecular_mass, material_unit_cell_volume, material_row["Z_param"])
            
            # Find percent of incident beam that is scattered (assume no absorption)
            # Find percent of incident beam that is absorbed (assume no scattering)
            # Find can mass in grams
            can_percent_scatter_ann, can_percent_absorb_ann, can_mass_ann = can_annulus(scatter_xs = material_scatter_xs, absorb_xs = material_absorb_xs, theory_density = material_theory_density, unit_cell_volume = material_unit_cell_volume, height = can_height_cm, inner_radius = can_R3_cm, outer_radius = new_can_R4_cm, volume = can_material_volume_cm3)
            
            # Populate dictionaries with id as key
            can_percent_scatter[id] = round(can_percent_scatter_ann, 2)
            can_percent_absorb[id] = round(can_percent_absorb_ann, 2)
            can_mass[id] = round(can_mass_ann, 3)
        
    # Convert dictionaries to DataFrames
    # Create dictionary df_dict of DataFrames
    df_dict = {
        "sample_mass_g": pd.DataFrame(sample_mass.items(), columns = ["id", "sample_mass_g"]),
        "can_volume_cm3": pd.DataFrame(can_volume_dict.items(), columns = ["id", "can_volume_cm3"]),
        "sample_moles": pd.DataFrame(sample_moles.items(), columns = ["id", "sample_moles"]),
        "percent_scatter": pd.DataFrame(percent_scatter.items(), columns = ["id", "percent_scatter"]),
        "percent_absorb": pd.DataFrame(percent_absorb.items(), columns = ["id", "percent_absorb"]),
        "can_percent_scatter": pd.DataFrame(can_percent_scatter.items(), columns = ["id", "can_percent_scatter"]),
        "can_percent_absorb": pd.DataFrame(can_percent_absorb.items(), columns = ["id", "can_percent_absorb"]),
        "can_mass_g": pd.DataFrame(can_mass.items(), columns = ["id", "can_mass_g"])
    }

    # Extract first DataFrame from dictionary
    # Extract DataFrame sample_mass from dictionary
    df = df_dict["sample_mass_g"]

    # Merge DataFrames in dicitonary with left join on id
    for key in list(df_dict.keys())[1:]:
        df = df.merge(df_dict[key], on = "id", how = "left")

    # Merge with left join on id with cans_desc DataFrame
    df_concat = pd.merge(cans_desc, df, on = "id", how = "left")
    
    # Set conditions and corresponding flags
    conditions = [
        df_concat["percent_scatter"] > 10,
        df_concat["percent_absorb"] > 10,
        (df_concat["percent_absorb"] + df_concat["percent_scatter"]) > 10,
        df_concat["can_percent_scatter"] > df_concat["percent_scatter"]
    ]
    flags = [f"% scatter > 10", f"% absorb > 10", f"% scatter + % absorb > 10", f"% scatter can > % scatter"]

    # Add empty column
    df_concat["flag"] = ""

    # Fill column with flag based on conditions
    for condition, flag in zip(conditions, flags):
        df_concat["flag"] = np.where(condition, df_concat["flag"] + ", " + flag, df_concat["flag"])

    # Remove leading comma
    df_concat["flag"] = df_concat["flag"].str.lstrip(', ')

    # Replace empty strings with NaN if no conditions are true
    df_concat["flag"] = df_concat["flag"].replace("", np.nan)

    # Report sample moles as sample millimoles 
    df_concat["sample_moles"] = round(df_concat["sample_moles"] * 1000, 3)
    df_concat.rename(columns = {"sample_moles": "sample_mmoles"}, inplace = True)

    # Drop drawing_number column
    df_concat.drop("drawing_number", axis = 1, inplace = True)

    sample_txt = f"""
    =============================
        Sample Information
    =============================
    Formula Unit (fu):  {mantid_formula}
    Z:                {Z_param}
    Lattice constant (a): {a} Å
    Lattice constant (b): {b} Å
    Lattice constant (c): {c} Å
    Lattice constant (alpha): {alpha}°
    Lattice constant (beta): {beta}°
    Lattice constant (gamma): {gamma}°
    Packing fraction: {pack_fraction}
    Incident neutron energy:   {neutron_energy} meV
    Incident neutron wavelength: {round(neutron_wavelength, 3)} Å
    ============================
    """

    sample_can_txt = f"""
    =========================================
        Sample Can Independent Values
    =========================================
    Calculated unit cell volume (Å^3): {round(unit_cell_volume, 3)}
    Molecular mass (g/mol/fu/Z): {round(molecular_mass / float(Z_param), 3)}
    Molecular mass per unit cell (g/mol/unit cell): {round(molecular_mass, 3)}
    Incoherent scattering cross section (bn/fu): {round(incoh_scatter_xs, 3)}
    Coherent scattering cross section (bn/fu): {round(coh_scatter_xs, 3)}
    Total scattering cross section (bn/fu): {round(scatter_xs, 3)}
    Absorption cross section (bn/fu): {round(absorb_xs, 3)}
    Scattering penetration depth (cm): {round(scatter_depth, 3)}
    Absorption penetration depth (cm): {round(absorb_depth, 3)}
    Total peneteration depth (cm):  {round(total_depth, 3)}
    Theoretical sample density (g/cc): {round(theory_density, 3)}
    Vanadium equivelent incoherent scatterer (g_V/g_fu): {round(van_equiv_incoherent, 4)}
    ==========================================
    """

    # flag_txt = f"""
    # ==========================
    #     Warning Flags
    # ==========================
    # (*)    % Scatter of sample > 10

    # (**)   % Absorb. of sample > 10

    # (***)  (% Scatter of sample + % Absorb. of sample) > 10

    # (****) % Scatter of can > % Scatter of sample
    # ==========================
    # """

    description_txt = f"""
    ====================================================
                    Variable Description
    ====================================================
    Id, sample can unique identifier.  This label contains the material of the can (Al, Cu, or V).

    Description, description of the sample can geometry and dimensions.

    Sample_mass_g, approximate mass of sample in grams to fill the sample can for the given packing fraction.

    Can_volume_cm3, available can volume in cubic centimeters.

    Sample_mmoles, millimoles of sample formula units corresponding to the sample_mass_g column.

    Percent_scatter, percent of neutrons scattered by the sample based upon the total neutron scattering cross-section and the geometry of the sample can.

    Percent_absorb, percent of neutrons absorbed by the sample based upon the wavelength dependent neutron absorption and the geometry of the sample can.  Calculation does not account for any neutron absorption resonances.

    Can_percent_scatter, percent of neutrons scattered by only the empty sample can based upon the total neutron scattering cross-section of the material and the geometry of the sample can.

    Can_percent_absorb, percent of neutrons absorbed by only the empty sample can based upon the wavelength dependent neutron absorption and the geometry of the sample can.  Calculation does not account for any neutron absorption resonances.

    Flag, warning flag based upon percent scattered and percent absorbed calculations.
    ====================================================
    """

    # Obtain date and time
    current_date = datetime.datetime.now().strftime("%d-%m-%y-%H-%M-%S")
    
    # Set filenames for text and CSV returns
    txt_filename = f"{mantid_formula}_{current_date}.txt"
    csv_filename = f"{mantid_formula}_{current_date}.csv"

    try:
        with open(txt_filename, "w", encoding = "utf-8") as txt_file:
            txt_file.write(sample_txt)
            txt_file.write(neutron_txt)
            txt_file.write(sample_can_txt)
            # txt_file.write(flag_txt)
            txt_file.write(description_txt)
    except Exception as e:
        print(f"An error occurred while writing to text file: {e}")

    # Write output_df as csv to /csv_filename.csv
    output_df = df_concat.dropna(subset = ["sample_mass_g"])
    output_df.to_csv(csv_filename, index = False)