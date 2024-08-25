# add_can_dict.py

import re
import pandas as pd
import numpy as np

# Define add_can_dict() function
def add_can_dict(file_path, can_type, can_material, sample_height_mm, sample_width_mm = None, sample_thick_mm = None, can_thick_f1_mm = None, can_thick_f2_mm = None, sample_inner_radius_mm = None, can_outer_radius_mm = None, sample_outer_radius_mm = None, insert_inner_radius_mm = None):
    """
    Parameters:
    can_type: "flatPlate", "cyl", or "annulus"
    can_material: e.g., "Cu", "Al", "V"
    sample_height_mm: Height of sample in mm
    """
    # Define replace_decimal_with_p() function
    def replace_decimal_with_p(value):
        value_str = str(value) 
        new_value = re.sub(r'\.', "p", value_str)
        return new_value
    
    # If `can_type` parameter is flatPlate
    if can_type == "flatPlate":
        """
        Parameters:
        sample_width_mm: Width of sample in mm
        can_thick_f1_mm: Thickness of first face of can in mm
        can_thick_f2_mm: Thickness of second face of can in mm
        sample_thick_mm: Thickness of sample in mm
        """  
        
        # Set filepath and read CSV as DataFrame   
        df = pd.read_csv(file_path)
        
        # Find can width, total thickness, and height
        can_width_mm = sample_width_mm # in mm
        can_total_thick_mm = can_thick_f1_mm + can_thick_f2_mm # in mm
        can_height_mm = sample_height_mm # in mm
        
        # Find sample area and can area
        sample_area_mm2 = sample_height_mm * sample_width_mm # in mm^2
        can_area_mm2 = sample_area_mm2 # in mm^2
        
        # Find sample volume and can volume
        sample_volume_mm3 = sample_area_mm2 * sample_thick_mm # in mm^3
        can_volume_mm3 = sample_volume_mm3 # in mm^3

        # Find volume of can material
        can_material_volume_cm3 = (sample_area_mm2 * can_total_thick_mm) / 1000 # in cm^3

        # Set unique can ID
        unique = f"{replace_decimal_with_p(sample_thick_mm)}mm"
        id = f"{can_type}_{unique}_{can_material}"

        # Create row for new can
        new_row = {
            "id": id,
            "material": can_material,
            "sample_thick_mm": sample_thick_mm,
            "sample_height_mm": sample_height_mm,
            "sample_width_mm": sample_width_mm,
            "can_width_mm": can_width_mm,
            "can_thick_f1_mm": can_thick_f1_mm,
            "can_thick_f2_mm": can_thick_f2_mm,
            "can_total_thick_mm": can_total_thick_mm,
            "can_height_mm": can_height_mm,
            "sample_area_mm2": sample_area_mm2,
            "can_area_mm2": can_area_mm2,
            "sample_volume_mm3": sample_volume_mm3,
            "can_volume_mm3": can_volume_mm3,
            "can_material_volume_cm3": can_material_volume_cm3
        }

    # If `can_type` parameter is cylindrical
    elif can_type == "cyl":
        """
        Parameters:
        sample_inner_radius_mm: Inner radius of sample in mm
        can_outer_radius_mm: Outer radius of can in mm
        """
        # Set filepath and read CSV as DataFrame
        df = pd.read_csv(file_path)
        
        # Find can inner radius and height
        can_inner_radius_mm = sample_inner_radius_mm # in mm
        can_height_mm = sample_height_mm # in mm
        
        # Find sample volume and can volume
        sample_volume_mm3 = np.pi * sample_height_mm * (sample_inner_radius_mm**2) # in mm^3
        can_volume_mm3 = sample_volume_mm3 # in mm^3

        # Find can thickness
        can_thick_mm = can_outer_radius_mm - can_inner_radius_mm # in mm

        # Find volume of can material
        can_material_volume_cm3 = np.pi * can_height_mm * (can_outer_radius_mm**2 - can_inner_radius_mm**2) / 1000 # in cm^3
        
        # Find can diammeter 
        can_diam_mm = 2 * sample_inner_radius_mm # in mm

        # Set unique can ID
        unique = f"{replace_decimal_with_p(can_diam_mm)}mm_{replace_decimal_with_p(sample_height_mm)}mm"
        id = f"{can_type}_{unique}_{can_material}"

        # Create row for new can
        new_row = {
            "id": id,
            "material": can_material,
            "sample_inner_radius_mm": sample_inner_radius_mm,
            "can_inner_radius_mm": can_inner_radius_mm,
            "sample_height_mm": sample_height_mm,
            "can_height_mm": can_height_mm,
            "sample_volume_mm3": sample_volume_mm3,
            "can_volume_mm3": can_volume_mm3,
            "can_outer_radius_mm": can_outer_radius_mm,
            "can_thick_mm": can_thick_mm,
            "can_material_volume_cm3": can_material_volume_cm3
        }

    # If `can_type` parameter is annulus
    elif can_type == "annulus":
        """
        Parameters:
        sample_outer_radius_mm: Outer radius of sample in mm
        sample_inner_radius_mm: Inner radius of sample in mm
        can_outer_radius_mm: Outer radius of cylindrical can in mm
        insert_inner_radius_mm: Inner radius of annular insert in mm
        """
        # Set filepath and read CSV as DataFrame
        df = pd.read_csv(file_path)
        
        # Find can inner radius and thickness
        can_inner_radius_mm = sample_outer_radius_mm # in mm
        can_thick_mm = can_outer_radius_mm - can_inner_radius_mm # in mm
        
        # Find sample thickness and volume
        sample_thick_mm = sample_outer_radius_mm - sample_inner_radius_mm # in mm
        sample_volume_mm3 = sample_height_mm * np.pi * (sample_outer_radius_mm**2 - sample_inner_radius_mm**2) # in mm^3

        # Find insert outer radius and thickness
        insert_outer_radius_mm = sample_inner_radius_mm # in mm
        insert_thick_mm = insert_outer_radius_mm - insert_inner_radius_mm # in mm

        # Convert insert/can radius and sample height from mm to cm
        can_R1_cm = insert_inner_radius_mm/10 # in cm
        can_R2_cm = insert_outer_radius_mm/10 # in cm
        can_R3_cm = can_inner_radius_mm/10 # in cm
        can_R4_cm = can_outer_radius_mm/10 # in cm
        can_height_cm = sample_height_mm/10 # in cm

        # Find inner and outer volume of can material
        can_material_volume_inner_cm3 = can_height_cm * np.pi * (can_R2_cm**2 - can_R1_cm**2)
        can_material_volume_outer_cm3 = can_height_cm * np.pi * (can_R4_cm**2 - can_R3_cm**2)

        # Find volume of can material
        can_material_volume_cm3 = can_material_volume_outer_cm3 + can_material_volume_inner_cm3 # in cm^3

        # Find difference in can outer radius and new can outer radius
        d = np.sqrt(can_R2_cm**2 - can_R1_cm**2 + can_R4_cm**2) - can_R4_cm # in cm
        new_can_R4_cm = can_R4_cm + d # in cm

        # Find new volume and total volume of can material
        new_can_material_volume_outer_cm3 = can_height_cm * np.pi * (new_can_R4_cm**2 - can_R3_cm**2) # in cm^3
        new_can_material_total_volume_cm3 = can_material_volume_inner_cm3 + new_can_material_volume_outer_cm3 # in cm^3
        
        # Set unique can ID
        unique = f"{replace_decimal_with_p(sample_height_mm)}mm_{replace_decimal_with_p(insert_inner_radius_mm)}mm"
        id = f"{can_type}_{unique}_{can_material}"

        # Create row for new can
        new_row = {
            "id": id,
            "material": can_material,
            "sample_height_mm": sample_height_mm,
            "sample_inner_radius_mm": sample_inner_radius_mm,
            "sample_outer_radius_mm": sample_outer_radius_mm,
            "sample_thick_mm": sample_thick_mm,
            "sample_volume_mm3": sample_volume_mm3,
            "can_inner_radius_mm": can_inner_radius_mm,
            "can_outer_radius_mm": can_outer_radius_mm,
            "can_thick_mm": can_thick_mm,
            "insert_inner_radius_mm": insert_inner_radius_mm,
            "insert_outer_radius_mm": insert_outer_radius_mm,
            "insert_thick_mm": insert_thick_mm,
            "can_R1_cm": can_R1_cm,
            "can_R2_cm": can_R2_cm,
            "can_R3_cm": can_R3_cm,
            "can_R4_cm": can_R4_cm,
            "can_height_cm": can_height_cm,
            "can_material_volume_cm3": can_material_volume_cm3, 
            "can_material_volume_inner_cm3": can_material_volume_inner_cm3,
            "can_material_volume_outer_cm3": can_material_volume_outer_cm3,
            "d": d,
            "new_can_R4_cm": new_can_R4_cm,
            "new_can_material_volume_outer_cm3": new_can_material_volume_outer_cm3,
            "new_can_material_total_volume_cm3": new_can_material_total_volume_cm3
        }
    else:
        raise ValueError("Invalid geometry type. Please enter 'flatPlate', 'cyl', or 'annulus'.")
    
    # Append new row to DataFrame and write back to CSV
    df = pd.concat([df, pd.DataFrame([new_row])], ignore_index = True)
    df.to_csv(file_path, index = False)