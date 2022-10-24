# Solar-PTES

A program for modelling Pumped Thermal Energy Storage (PTES) and hybrid solar-PTES using MATLAB and CoolProp.

## Capabilities

This program can conduct:
- Design point cycle calculations
- Off-design calculations (variable power output, ambient temperature, storage tank temperature)
- Economic estimates (capital cost, LCOS) and uncertainty analysis
- Parametric studies over two variables
- Multi-objective optimization

The majority of recent development has concentrated on PTES systems using liquid thermal storage (such as molten salts) and Brayton cycles (e.g. working fluids such as air, argon, nitrogen).

Previously, many other variations have been modelled within this program framework, including PTES using supercritical-CO2, steam-Rankine cycles, and packed bed thermal energy storage. Previous work also developed "solar-PTES" concepts where PTES was hybridized with Concentrating Solar Power (CSP).


## Simple Interface

The program has many capabilities and many of them are (currently) accessed in a non-intuitive way. 

The so-called Simple Interface provides an introduction to running the code and enables the user to vary some key parameters, such as turbo-machinery and heat exchanger performance and system temperatures. The parameters available to edit largely correspond to those that define the PTES system modelled in the System Advisor Model (beta version).

To run the program using the Simple Interface:

1. At the start of `PTES_code\PTES_scripts\INPUTS.m` ensure that `simple_interface = 1;`
2. Open `PTES_code\PTES_scripts\SET_SIMPLE_INTERFACE.m` and edit the variables as desired.
3. Set working directory to `PTES_code`
4. Type `main` in the MATLAB command window

Once the program has completed running, several useful outputs are provided. Information regarding the system performance is written to the command window. 

## Integration with SAM

PTES models are being integrated into the System Advisor Model (SAM) developed by NREL (<https://sam.nrel.gov/>). (Currently only available as a beta release). 

Key information from the Simple Interface MATLAB model is written into a .json file which can be found in `PTES_code\Outputs\PTES_output_for_SAM.json`. This json can be read into SAM using a macro. 

## Accessing more advanced capabilities

Most of the capabilities mentioned above can only be accessed if the Simple Interface is turned off (`simple_interface = 0;`).

The user can then specify many options using the input files `PTES_code\PTES_scripts\INPUTS.m` and `PTES_code\PTES_scripts\JB_RANK_INPUTS.m`

# References

1. J. D. McTigue, P. Farres-Antunez, K. Sundarnath, C. N. Markides, A. J. White, "Techno-economic analysis of recuperated Joule-Brayton pumped thermal energy storage systems", Energy Conversion & Management, vol. 252, 115016, 2022, <https://doi.org/10.1016/j.enconman.2021.115016>

2. J. D. McTigue, P. Farres-Antunez, C. N. Markides, A. J. White, "Pumped thermal energy storage with liquid storage", chapter in Encyclopedia of Energy Storage, Elsevier, 2021, <https://doi.org/10.1016/B978-0-12-819723-3.00054-8>

3. J. Martinek, J. Jorgenson, J. D. McTigue, "On the operational characteristics and economic value of pumped thermal energy storage", Journal of Energy Storage, vol. 52, 105005, 2022, <https://doi.org/10.1016/j.est.2022.105005>

4. J. D. McTigue, P. Farres-Antunez, C. N. Markides, A. J. White, "Integration of heat pumps with solar thermal systems for energy storage", chapter in Encyclopedia of Energy Storage, Elsevier, 2021, <https://doi.org/10.1016/B978-0-12-819723-3.00067-6>

5. J. D. McTigue, P. Farres-Antunez, T. Neises, A. J. White, "Supercritical CO2 heat pumps and power cycles for concentrating solar power", SolarPACES conference, virtual conference, October 2020, <https://www.nrel.gov/docs/fy21osti/77955.pdf>

6. I. H. Bell, J. Wronski, S. Quoilin, V. Lemort, "Pure and Pseudo-pure Fluid Thermophysical Property Evaluation and the Open-Source Thermophysical Property Library CoolProp", Ind. Eng. Chem. Res. 2014, 53, 6, 2498–2508, DOI: 10.1021/ie4033999



