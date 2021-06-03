"""
/*---------------------------------------------------------------------------*\
 ####                 ####     #   #             |
 #                    #  #     ##  #             | FoamPython
 ##  #### ####  ##### #### # # #   #### #### ### | v1.0
 #   #  # #  #  # # # #     #  # # #  # #  # # # |
 #   #### ##### # # # #    #    #  #  # #### # # |
-------------------------------------------------------------------------------

Author
	Robert Anderluh, 2021
	
Description
    Include file for the cfdTools folder.
    
\*---------------------------------------------------------------------------*/
"""

import os

# Save current directory
oldDir = os.getcwd()

# Get the directory of this file
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# The location of continuityErrs code
continuityErrs = open("incompressible/continuityErrs.py").read()

# Go back to the old directory
os.chdir(oldDir)

# Import the constrainHbyA function
from src.finiteVolume.cfdTools.general.constrainHbyA.constrainHbyA import *

"""
// ************************************************************************* //
"""
