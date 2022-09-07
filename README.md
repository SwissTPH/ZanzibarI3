# Zanzibar malaria transmission model for imported, introduced and indigenous cases

Copyright (C) 2022 Swiss Tropical and Public Health Institute

This malaria model is free software; you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.

This model is built for studying imported, introduced and indigenous cases in settings approaching malaria elimination. It includes parameters for studying reactive case detection, treatment or infection prevention in travellers, and general reductions in transmission rates.

### Data
* 'RADZEC_data.csv': point value estimates of relevant parameter values for running the model
* 'distributions.csv': RADZEC data values used to parameterise distributions for the uncertainty analysis
* 'make_param_sets.py': sample from the distributions representing the uncertainty in the RADZEC data for the uncertainty analysis
* 'param_100.csv': parameter values drawn and used to conduct the uncertainty analysis

### Interventions
Parameter values for the various intervention scenarios can be found here.
A script for making the intervention tables is also included ('make_interventions.py').

### Code
All code to run the model in Python can be found in the 'code' folder. Separate files are provided for different interventions and analyses.

### Plotting
Code for postprocessing and plotting can be found here.

### Sensitivity Analysis
The code for samplying for the sensitivity analysis, running the model, and then calculating the Sobol indices from the inputs and outputs is included here.
