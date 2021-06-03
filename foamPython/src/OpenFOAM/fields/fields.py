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
    Generic fields class. The foundation of the fields class hierarchy
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np
import os
import shutil

from src.OpenFOAM.fields.writeHeader import HEADER


class fields():
	
	fieldTypeName_ = "fieldType"
	className_ = "fields"
	noComponents_ = 0
	componentsNames_ = None
	
	"""------------------------ General functions ---------------------------"""
	
	def updateBoundaryDefinition(self, boundaryDefinition):
		
		self.data.boundaryDefinition_ = boundaryDefinition
		
	def updateFieldName(self, fieldName):
		
		self.data.fieldName_ = fieldName
	
"""
// ************************************************************************* //
"""
