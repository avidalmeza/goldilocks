# Goldilockσ: The "Just Right" Can
Goldilockσ calculates the total scattering and absorption cross-sections for all standard powder cans (flat plate, cylinder, and annular) available at the Spallation Neutron Source and High Flux Isotope Reactor.

## Parameters
### Sample
Users can set a sample with a Crystallographic Information File (CIF) or molecular formula.
#### Crystallographic Information File
A CIF is a standard text file format for representing crystallographic information. If your sample does not include an isotope, you can enter your CIF as is. 

If your sample includes an isotope, open your CIF in a text editor and modify it to match the accepted format. Users need to edit the `_chemical_formula_sum` variable in the CIF to vary composition or include isotopes.

For example:
> _chemical_formula_sum 'Cl1.45 Li5.47 P1 S4.55'

> _chemical_formula_sum 'Cl1.45 (Li7)5.47 P1 S4.55'

For example: 
> _chemical_formula_sum 'Ba6 Co3 La6 O36 Te6 H3'

> _chemical_formula_sum 'Ba6 Co3 La6 O36 Te6 H0.1 (H2)2.9'

#### Molecular Formula
Refer to the [Mantid project's documentation](https://docs.mantidproject.org/v6.6.0/concepts/Materials.html) for the accepted format of a molecular formula. Each element in a molecular formula is followed by the number of atoms for that element, specified without a hyphen because each element is separated from other elements using a hyphen. Isotopes must be enclosed by parenthesis.

For example:
>'Ba6-Co3-La6-O36-Te6-H0.1-(H2)2.9'

### Lattice Parameters
Enter the values of the lattice parameters a, b, c, alpha, beta, gamma, where a, b, c are lengths of the sides of the unit cell and alpha is the angle between b and c, beta is the angle between a and c, and gamma is the angle between a and b. Enter the lattice constants (a, b, c) in Angstroms. Enter the lattice angles (alpha, beta, and gamma) in degrees. Users who enter a CIF do not need to enter lattice parameters.
### Formula Units
Enter the number of formula units (Z) in a crystallographic unit cell to calculate unit cell volume. Z is a unitless parameter. Users who enter a CIF or enter a molecular formula for a new molecule or compound do not need to enter this parameter.
### Incident Neutron Energy
Enter the incident (incoming) neutron energy in meV. Users must set value between 0.01 to 100000 meV. Default value is 14.7 meV. Users can enter either the incident neutron energy or neutron wavelength.
### Neutron Wavelength
Enter the neutron wavelength in Angstroms. Default value is 2.359 Å. Users can enter either the incident neutron energy or neutron wavelength.
### Packing Fraction
Enter the fraction of volume in the structure that is occupied by constituent particles. It is the ratio of the effective density of the material in the container to the calculated density of the material if it were a solid. The packing fraction is dimensionless and the value should be less than or equal to one. The packing fraction for typical powders ranges from 0.3 to 0.7. Default value is set to 1.0 and assumes a solid sample.
### Density
Enter the sample density in units of grams per cubic centimeter. Only users who enter a molecular formula for a new molecule or compound need to enter this parameter.

## Returns
- **sample_mass_g** (float): Approximate mass of sample in grams to fill the sample can for the given packing fraction
- **can_volume_cm3** (float): Available can volume in cubic centimeters (cm^3)
- **sample_mmoles** (float): Millimoles of sample formula units corresponding to approximate mass of sample
- **percent_scatter** (float): Percent of neutrons scattered by the sample based upon the total neutron scattering cross-section and the geometry of the sample can
- **percent_absorb** (float): Percent of neutrons absorbed by the sample based upon the wavelength dependent neutron absorption and the geometry of the sample can; does not account for any neutron absorption resonances
- **can_percent_scatter** (float): Percent of neutrons scattered by only the empty sample can based upon the total neutron scattering cross-section of the material and the geometry of the sample can
- **can_percent_absorb** (float): Percent of neutrons absorbed by only the empty sample can based upon the wavelength dependent neutron absorption and the geometry of the sample can; does not account for any neutron absorption resonances
- **flag** (str): Warning flag based upon percent scattered and percent absorbed calculations

## Updates
Goldilockσ is built to maintain scientific reproducibility and leverages local dictionaries and the Mantid project's Python API.

Goldilockσ provides estimates for all standard powder cans (flat plate, cylinder, and annular) available at SNS and HFIR. If there is a new can of these geometries, users can edit the corresponding CSV file ([flatPlate.csv](/src/dict/flatPlate.csv), [cylindrical.csv](/src/dict/cylindrical.csv), and [annulus.csv](/src/dict/annulus.csv)) with [add_can.ipynb](/src/add_can.ipynb). Alternatively, users can manually update the corresponding CSV file. In both cases, users must manually add an entry for the can description to [canDescription.csv](/src/dict/canDescription.csv).

Goldilockσ provides estimates for powder cans made of aluminum, copper, or vanadium. In the case of a new material, users can edit [material.csv](/src/dict/material.csv) through [add_can.ipynb](/src/add_can.ipynb) with a CIF. Alternatively, users without a CIF can manually set the values of the lattice parameters and formula units in [material.csv](/src/dict/material.csv).