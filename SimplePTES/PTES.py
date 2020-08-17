# This is a python program to provide a simple model of Pumped Thermal Energy Storage (PTES)

# The PTES is a recuperated Joule-Brayton cycle with liquid storage tanks.
# See McTigue 2015, Laughlin 2019 for more details.

# Code developed by Josh McTigue, July 2020 for Jen King's NREL LDRD
# JoshuaDominic.McTigue@nrel.gov

import PTES_functions as PT
import importlib as ilib
ilib.reload(PT) # This allows the module to be called repeated from the interpreter even if changes are being made to the module

# Load the inputs
inp = PT.inputs()
print('Loaded inputs')

# Set up PTES cycle class
CYC = PT.Cycle(inp)

# Call design charging cycle
CYC = PT.charge(inp,CYC)

# Call design discharging cycle
CYC = PT.discharge(inp,CYC)

# Calculate design point performance
CYC = PT.Cycle.performance(CYC)

