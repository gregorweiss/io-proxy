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
    Surface vector field class
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.surfaceField.surfaceField import *


class surfaceVectorField(surfaceField):
	
	fieldTypeName_ = "vector"
	className_ = "surfaceVectorField"
	noComponents_ = 3
	componentsNames_ = ("x", "y", "z")
	

"""
// ************************************************************************* //
"""
