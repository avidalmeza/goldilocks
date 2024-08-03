# add_can_dict.py

import os
import re
import numpy as np
import pandas as pd
import csv

# Obtain current working directory filepath
thisDir = os.getcwd()

# Read sample can flatPlate csv
flatPlate_df = pd.read_csv(os.path.join(thisDir, 'src', 'dict', 'flatPlate.csv'))
# Read sample can cylindrical csv
cyl_df = pd.read_csv(os.path.join(thisDir, 'src', 'dict', 'cylindrical.csv'))
# Read sample can annulus csv
annulus_df = pd.read_csv(os.path.join(thisDir, 'src', 'dict', 'annulus.csv'))

# Define replace_decimal_with_p() function
def replace_decimal_with_p(value):
    # Convert float to string
    value_str = str(value)
    # Replace decimal point with 'p'
    new_value = re.sub(r'\.', 'p', value_str)
    return new_value

def add_can_dict(can_type, can_material, sample_height_mm, sample_width_mm = None, sample_thick_mm = None, 
                 can_thick_f1_mm = None, can_thick_f2_mm = None, sample_inner_radius_mm = None, can_outer_radius_mm = None, 
                 sample_outer_radius_mm = None, insert_inner_radius_mm = None):
    
    if can_type == 'flatPlate':
        # sample_thick_mm, sample_height_mm, sample_width_mm, can_thick_f1_mm, can_thick_f2_mm

        file_path = os.path.join(thisDir, 'src', 'dict', 'flatPlate.csv')
        
        can_width_mm = sample_width_mm
        can_total_thick_mm = can_thick_f1_mm + can_thick_f2_mm
        can_height_mm = sample_height_mm
        
        sample_area_mm2 = sample_height_mm*sample_width_mm
        can_area_mm2 = sample_area_mm2
        
        sample_volume_mm3 = sample_area_mm2*sample_thick_mm
        can_volume_mm3 = sample_volume_mm3

        can_material_volume_cm3 = (sample_area_mm2*can_total_thick_mm)/1000

        unique = f'{replace_decimal_with_p(sample_thick_mm)}'
        id = f'{can_type}_{unique}_{can_material}'

        new_row = {'id': id,
                   'material': can_material,
                   'sample_thick_mm': sample_thick_mm,
                   'sample_height_mm': sample_height_mm,
                   'sample_width_mm': sample_width_mm,
                   'can_width_mm': can_width_mm,
                   'can_thick_f1_mm': can_thick_f1_mm,
                   'can_thick_f2_mm': can_thick_f2_mm,
                   'can_total_thick_mm': can_total_thick_mm,
                   'can_height_mm': can_height_mm,
                   'sample_area_mm2': sample_area_mm2,
                   'can_area_mm2': can_area_mm2,
                   'sample_volume_mm3': sample_volume_mm3,
                   'can_volume_mm3': can_volume_mm3,
                   'can_material_volume_cm3': can_material_volume_cm3
        }

        # Insert `new_row` DataFrame to end of `flatPlate_df` DataFrame and reset index
        flatPlate_df.loc[len(flatPlate_df)] = new_row
        flatPlate_df = flatPlate_df.reset_index(drop = True)

        # Write `flatPlate_df` as csv to src/dict/flatPlate.csv
        flatPlate_csv_filepath = os.path.join(thisDir, 'src', 'dict', 'flatPlate.csv')
        flatPlate_df.to_csv(flatPlate_csv_filepath, index = False)

    elif can_type == 'cyl':
        # sample_inner_radius_mm, sample_height_mm, can_outer_radius_mm

        file_path = os.path.join(thisDir, 'src', 'dict', 'cylindrical.csv')
        
        can_inner_radius_mm = sample_inner_radius_mm
        can_height_mm = sample_height_mm
        
        sample_volume_mm3 = np.pi*sample_height_mm*(sample_inner_radius_mm**2)
        can_volume_mm3 = sample_volume_mm3

        can_thick_mm = can_outer_radius_mm - can_inner_radius_mm

        can_material_volume_cm3 =  np.pi*can_height_mm*(can_outer_radius_mm**2 - can_inner_radius_mm**2)/1000

        unique = f'{replace_decimal_with_p(sample_inner_radius_mm)}_{replace_decimal_with_p(sample_height_mm)}'
        id = f'{can_type}_{unique}_{can_material}'

        # Create `new_row` DataFrame 
        new_row = {'id': id,
                   'material': can_material,
                   'sample_inner_radius_mm': sample_inner_radius_mm,
                   'can_inner_radius_mm': can_inner_radius_mm,
                   'sample_height_mm': sample_height_mm,
                   'can_height_mm': can_height_mm,
                   'sample_volume_mm3': sample_volume_mm3,
                   'can_volume_mm3': can_volume_mm3,
                   'can_outer_radius_mm': can_outer_radius_mm,
                   'can_thick_mm': can_thick_mm,
                   'can_material_volume_cm3': can_material_volume_cm3
        }
        
        # Insert `new_row` DataFrame to end of `cyl_df` DataFrame and reset index
        cyl_df.loc[len(cyl_df)] = new_row
        cyl_df = cyl_df.reset_index(drop = True)

        # Write `cyl_df` as csv to src/dict/cylindrical.csv
        cyl_csv_filepath = os.path.join(thisDir, 'src', 'dict', 'cylindrical.csv')
        cyl_df.to_csv(cyl_csv_filepath, index = False)

    elif can_type == 'annulus':
        # sample_height_mm, sample_inner_radius_mm, sample_outer_radius_mm, can_outer_radius_mm, insert_inner_radius_mm

        file_path = os.path.join(thisDir, 'src', 'dict', 'annulus.csv')
        
        sample_thick_mm = sample_outer_radius_mm - sample_inner_radius_mm
        
        sample_volume_mm3 = sample_height_mm*np.pi*(sample_outer_radius_mm**2 - can_inner_radius_mm**2)/1000
        
        can_inner_radius_mm = sample_outer_radius_mm
        can_thick_mm = can_outer_radius_mm - can_inner_radius_mm
        
        insert_outer_radius_mm = sample_inner_radius_mm
        insert_thick_mm = insert_outer_radius_mm - insert_inner_radius_mm

        can_R1_cm = insert_inner_radius_mm/10
        can_R2_cm = insert_outer_radius_mm/10
        can_R3_cm = can_inner_radius_mm/10
        can_R4_cm = can_outer_radius_mm/10
        can_height_cm = sample_height_mm/10

        can_material_volume_inner_cm3 = can_height_cm*np.pi*(can_R2_cm**2 - can_R1_cm**2)
        can_material_volume_outer_cm3 = can_height_cm*np.pi*(can_R4_cm**2 - can_R3_cm**2)

        can_material_volume_cm3 = (can_material_volume_outer_cm3 + can_material_volume_inner_cm3)/1000

        d = np.sqrt(can_R2_cm**2 - can_R1_cm**2 + can_R4_cm**2) - can_R4_cm # in cm
        new_can_R4_cm = can_R4_cm + d # in cm

        new_can_material_volume_outer_cm3 = can_height_cm*np.pi*(new_can_R4_cm**2 - can_R3_cm**2) # in cm3
        new_can_material_total_volume_cm3 = can_material_volume_inner_cm3 + new_can_material_volume_outer_cm3 # in cm3

        unique = f'{replace_decimal_with_p(sample_height_mm)}_{replace_decimal_with_p(insert_inner_radius_mm)}'
        id = f'{can_type}_{unique}_{can_material}'

        # Create `new_row` DataFrame 
        new_row = {'id': id,
                   'material': can_material,
                   'sample_height_mm': sample_height_mm,
                   'sample_inner_radius_mm': sample_inner_radius_mm,
                   'sample_outer_radius_mm': sample_outer_radius_mm,
                   'sample_thick_mm': sample_thick_mm,
                   'sample_volume_mm3': sample_volume_mm3,
                   'can_inner_radius_mm': can_inner_radius_mm,
                   'can_outer_radius_mm': can_outer_radius_mm,
                   'can_thick_mm': can_thick_mm,
                   'insert_inner_radius_mm': insert_inner_radius_mm,
                   'insert_outer_radius_mm': insert_outer_radius_mm,
                   'insert_thick_mm': insert_thick_mm,
                   'can_R1_cm': can_R1_cm,
                   'can_R2_cm': can_R2_cm,
                   'can_R3_cm': can_R3_cm,
                   'can_R4_cm': can_R4_cm,
                   'can_height_cm': can_height_cm,
                   'can_material_volume_cm3': can_material_volume_cm3, 
                   'can_material_volume_inner_cm3': can_material_volume_inner_cm3,
                   'can_material_volume_outer_cm3': can_material_volume_outer_cm3,
                   'd': d,
                   'new_can_R4_cm': new_can_R4_cm,
                   'new_can_material_volume_outer_cm3': new_can_material_volume_outer_cm3,
                   'new_can_material_total_volume_cm3': new_can_material_total_volume_cm3
        }

        # Insert `new_row` DataFrame to end of `annulus_df` DataFrame and reset index
        annulus_df.loc[len(annulus_df)] = new_row
        annulus_df = annulus_df.reset_index(drop = True)

        # Write `annulus_df` as csv to src/dict/annulus.csv
        annulus_csv_filepath = os.path.join(thisDir, 'src', 'dict', 'annulus.csv')
        annulus_df.to_csv(annulus_csv_filepath, index = False)

    else:
        print("Invalid geometry type. Please enter 'flatPlate', 'cyl', or 'annulus'.")
        return