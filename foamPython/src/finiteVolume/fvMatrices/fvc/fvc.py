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
    Finite volume matrix contributions, explicit
    
    Note: fvc in OpenFOAM returns a field, as opposed to this implementation,
    where fvc returns an fvMatrix
    
\*---------------------------------------------------------------------------*/
"""

from src.finiteVolume.fvMatrices.fvc.fvSchemes.fvSchemes import *
from src.finiteVolume.fvMatrices.fvMatrix import *

class fvc(fvMatrix, fvSchemes):

	"""------------- Matrix components definition functions ----------------"""
	
	@classmethod
	def defineMatrix(self, psi, operator, scheme, fvVariables):
		
		mesh = psi.data.mesh_
		
		# The dictionary of operator classes
		operatorsDict = self.operatorsDict_
		
		# Assign the appropriate operator class
		try:
			operatorClass = operatorsDict[operator]
		except:
			raise RuntimeError("Operator " + operator + " is not implemented!")
		
		# Assign the appropriate function inside the operator class
		operatorSchemeFunc = eval("operatorClass." + scheme)
		
		# Call the appropriate function
		source, lower, diag, upper = operatorSchemeFunc(mesh, psi, fvVariables)
		
		return source, lower, diag, upper

"""
// ************************************************************************* //
"""
