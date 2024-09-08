# Goldilockσ: The “Just Right” Can 
<p align='center' width='100%'><img width='33%' src='/src/images/goldilocks.png'></p>

Neutrons are [a non-destructive probe to study materials on an atomic scale](https://neutrons.ornl.gov/industry/why-neutrons). [Because of their neutral charge](https://cen.acs.org/articles/88/i8/Making-Use-Neutrons.html), neutrons can penetrate deep into a sample [across a range of sample environments](https://www.isis.stfc.ac.uk/Pages/Why-and-how-to-use-neutrons-and-muons.aspx). Additionally, [due to their mass and nuclear spin](https://cen.acs.org/articles/88/i8/Making-Use-Neutrons.html), neutrons are also sensitive to the location and orientation of magnetic moments. Unlike other probes, neutrons scatter from the nuclei of atoms rather than the electron clouds, making them sensitive to light elements and allowing them to [distinguish between different isotopes](https://doi.org/10.2138/gselements.17.3.155). Neutrons' wavelength are also [sensitive to atomic length scales](https://www.ornl.gov/blog/what-makes-neutron-scattering-unique) and useful for studying the atomic structure of a material. As a result, neutron sciences research offers invaluable insight into challenges across energy, nanotechnology, transportation, communication, and several other areas. 

Oak Ridge National Laboratory (ORNL) operates 11 spectrometers and 12 diffractometers at the Spallation Neutron Source (SNS) and High Flux Isotope Reactor (HFIR). These instruments help determine the magnetic and crystallographic structure of materials and [reveal how to make the materials stronger, lighter, and perform better](https://neutrons.ornl.gov/industry/why-neutrons). Sample preparation, including the selection of the right sample holder, is key to a successful experiment. If a sample can is too small, then there is insufficient scattering and the measurement is too slow. In which case, there's not enough data. However, if a sample can is too large, surplus scattering occurs and the measurement may contain distortions and be inaccurate due to multiple scattering events of the neutron. Likewise, if there is too little scattering, the experiment will require too much time. At very long wavelengths, for certain materials, the scattering from the sample environment and sample can will deplete the beam and make measurements more lengthy.

The cross-section ($\sigma$), measured in barns, is the likelihood of the incident neutron interacting with a target nucleus and differs across elements and isotopes as well as sample holder geometries. Goldilockσ calculates the total scattering and absorption cross-sections for all standard powder cans (flat plate, cylinder, and annular). It provides a systematic and reproducible way to decide which can is the best can for neutron scattering measurements. So that external users and internal scientists at SNS and HFIR can choose the "just right" can for their experiment.

## Instructions
Users can clone this GitHub repository to run the calculator.

1. Clone the project from the command line
```
git clone https://github.com/avidalmeza/goldilocks.git
```

2. Create an environment from [*environment.yml*](/environment.yml) and activate
```
conda env create -f environment.yml
conda activate goldilocks
```

Activating the virtual environment is necessary to avoid issues with dependencies.

### Pathways 
Users can run the calculator from within a Jupyter Notebook, the Terminal (Linux/macOS), Command Prompt (Windows), or any other Command Line Interface (CLI).

#### Jupyter Notebook
3. Open [*xs_calculator.ipynb*](/xs_calculator.ipynb) and set calculator inputs

4. Run the calculator directly within the Jupyter Notebook

#### Terminal/Command Prompt
3. Open [*xs_calculator.py*](/xs_calculator.py) and set calculator inputs

4. Use CLI of choice to execute calculator
```
python xs_calculator.py
```

Refer to the [documentation](/documentation.md) for detailed instructions on calculator inputs and usage guidelines.

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

## Acknowledgments
This project was supported in part by an appointment to the ORNL GEM Fellow Internship Program, sponsored by the U.S. Department of Energy and administered by the Oak Ridge Institute for Science and Education.

### Authors and Contributors 
- Alessandra Vidal Meza, University of California, Santa Barbara
- Matthew Stone, Neutron Scattering Division, ORNL
- Andrew Christianson, Materials Science and Technology Division, ORNL
- Yuanpeng Zhang, Neutron Scattering Division, ORNL

### Artwork Credits
- *Neutron Scattering Barns* by Milo Stone