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
    Class which contains the temporal contributions for fvm
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

class ddtClass:
	
	class ddt:
		
		def Euler(mesh, psi, fvVariables):
			
			deltaT = fvVariables[0]
			cellVolumes = mesh.geometryData.V_
			cellValues = psi.data.cellValues_
			
			noComponents = psi.noComponents_
			
			# Total number of cells
			nCells = mesh.read("geometryData.meshSize")
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			source = np.zeros((noComponents, nCells), dtype = float)
			
			lower = np.zeros((noComponents, nInternalFaces), dtype = float)
			diag = np.zeros((noComponents, nCells), dtype = float)
			upper = np.zeros((noComponents, nInternalFaces), dtype = float)
			
			"""SOURCE DEFINITION"""
			# b = V * psi_old / deltaT
			for cmpt in range(noComponents):
				for cellIndex in range(nCells):
					source[cmpt][cellIndex] = 								   \
						cellVolumes[cellIndex] * cellValues[cmpt][cellIndex]   \
						/ deltaT
						
			"""A MATRIX DEFINITION"""
			# ap = V / deltaT
			for cmpt in range(noComponents):
				for cellIndex in range(nCells):
					diag[cmpt][cellIndex] = 								   \
						cellVolumes[cellIndex] 								   \
						/ deltaT
			
			return source, lower, diag, upper

"""
// ************************************************************************* //
"""
