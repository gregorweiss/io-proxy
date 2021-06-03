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
    Class which contains the functions needed to assemble the explicit fvMatrix,
    fvc
    
\*---------------------------------------------------------------------------*/
"""

from src.finiteVolume.fvMatrices.fvc.fvSchemes.grad.grad import *
from src.finiteVolume.fvMatrices.fvc.fvSchemes.div.div import *

class fvSchemes():
	
	# Dictionary of implemented operators
	operatorsDict_ = {														   \
					  "grad": gradClass.grad,								   \
					  "div": divClass.div									   \
					 }

"""
// ************************************************************************* //
"""
