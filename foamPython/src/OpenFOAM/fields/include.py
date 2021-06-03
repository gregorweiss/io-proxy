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
    Include file for the fields class hierarchy
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.volField.volScalarField.volScalarField import *
from src.OpenFOAM.fields.volField.volVectorField.volVectorField import *

from src.OpenFOAM.fields.surfaceField.surfaceScalarField.surfaceScalarField import *
from src.OpenFOAM.fields.surfaceField.surfaceVectorField.surfaceVectorField import *

"""
// ************************************************************************* //
"""
