# Goldilockσ: The "Just Right" Can
Goldilockσ calculates the total scattering cross-section and absorption cross-section for all standard powder cans (flat plate, cylinder, and annular) and how much scattering is due to the can itself. 

## Parameters
### Sample
Users can set a sample with a Crystallographic Information File (CIF) or molecular formula.
#### Crystallographic Information File
A CIF is a standard text file format for representing crystallographic information. If your sample does not include an isotope, you can enter your CIF as is. If your sample includes an isotope, open your CIF in a text editor and modify it to match the accepted format.
#### Molecular Formula
Refer to the [Mantid project's documentation](https://docs.mantidproject.org/v6.6.0/concepts/Materials.html) for the accepted format of a molecular formula. Each element in a molecular formula is followed by the number of the atoms for that element, specified without a hyphen, because each element is separated from other elements using a hyphen. Isotopes must be enclosed by parenthesis.

### Lattice Parameters
Enter the values of the lattice parameters a, b, c, alpha, beta, gamma, where a, b, c are lengths of the sides of the unit cell and alpha is the angle between b and c, beta is the angle between a and c, and gamma is the angle between a and b. Enter a, b, c in Angstroms. Enter alpha, beta, and gamma in degrees. Users who enter a CIF do not need to enter lattice parameters.

### Formula Units
Enter the number of formula units (Z) in a crystallographic unit cell to calculate unit cell volume. Z is a unitless parameter. Users who enter a CIF do not need to enter this parameter.

### Incident Neutron Energy
Enter the incident (incoming) neutron energy in meV. Users must set value between  0.01 to 100000 meV.

### Packing Fraction
Enter the fraction of volume in the structure that is occupied by constituent particles. The packing fraction is dimensionless.

## Returns
- **sample_mass_g** (float): Mass of the sample in grams
- **can_volume_cm3** (float): Volume of the can in cubic centimeters (cm^3)
- **sample_moles** (float): Number of moles of the formula unit in the sample
- **percent_scatter** (float): Percentage of the incident beam that is scattered, with no absorption, for the sample
- **percent_absorb** (float): Percentage of the incident beam that is absorbed, with no scattering, for the sample
- **can_percent_scatter** (float): Percentage of the incident beam that is scattered, with no absorption, for the can
- **can_percent_absorb** (float): Percentage of the incident beam that is absorbed, with no scattering, for the can
- **can_mass_g** (float): Mass of the can in grams

## Updates
Goldilockσ is built to mantain scientific reproducibility and leverages local dictionaries in addition to the Mantid project's Python API.
### New Sample Can
Goldilockσ provides estimates for all standard powder cans (flat plate, cylinder, and annular) available at SNS and HFIR. 

In the case of a new can, users can edit the corresponding CSV file ([flatPlate.csv](/src/dict/flatPlate.csv), [cylindrical.csv](/src/dict/cylindrical.csv), and [annulus.csv](/src/dict/annulus.csv)) with [add_can_dict.py](/src/add_can_dict.py). Alternatively, users can manually update the corresponding CSV file and manually add an entry for the can description to [canDescription.csv](/src/dict/canDescription.csv).
### New Sample Can Material
Goldilockσ provides estimates for powder cans made of aluminium, copper, or vanadium. In the case of a new sample can material, users can edit [material.csv](/src/dict/material.csv) through [add_material_dict.py](/src/add_material_dict.py) with a CIF. Alternatively, users without a CIF can set the values manually of the lattice parameters and formula units in [material.csv](/src/dict/material.csv).