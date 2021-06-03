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
    Class which contains the functions needed to assemble the implicit fvMatrix,
    fvm
    
\*---------------------------------------------------------------------------*/
"""

from src.finiteVolume.fvMatrices.fvm.fvSchemes.laplacian.laplacian import *
from src.finiteVolume.fvMatrices.fvm.fvSchemes.ddt.ddt import *
from src.finiteVolume.fvMatrices.fvm.fvSchemes.div.div import *

class fvSchemes(laplacianClass, ddtClass, divClass):
	
	# Dictionary of implemented operators
	operatorsDict_ = {														   \
					  "laplacian": laplacianClass.laplacian,				   \
					  "ddt" : ddtClass.ddt,									   \
					  "div" : divClass.div									   \
					 }

"""
// ************************************************************************* //
"""
