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
    Class which contains the divergence contributions for fvm
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

from src.finiteVolume.fvMatrices.fvm.fvSchemes.div.boundaryContributions import *

class divClass:
	
	class div(divBoundaryContributions):
		
		@classmethod
		def upwind(self, mesh, psi, fvVariables):
			
			phi = fvVariables[0].data.faceValues_
			Sf = mesh.geometryData.Sf_
			C = mesh.geometryData.C_
			
			noComponents = psi.noComponents_
			
			# Total number of cells
			nCells = mesh.read("geometryData.meshSize")
			# Total number of faces
			nFacesTot = np.size(mesh.connectivityData.owner_)
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			source = np.zeros((noComponents, nCells), dtype = float)
			
			lower = np.zeros((noComponents, nInternalFaces), dtype = float)
			diag = np.zeros((noComponents, nCells), dtype = float)
			upper = np.zeros((noComponents, nInternalFaces), dtype = float)
					 
					 
			"""A MATRIX DEFINITION"""
			# Add internal faces contributions
			for cmpt in range(noComponents):
				
				for faceIndex in range(nInternalFaces):
				
					# aP_owner = max(F, 0)
					# aN_owner = min(F, 0)
					# aP_neighbour = max(-F, 0)
					# aN_neighbour = min(-F, 0)
					
					cellIndexP = mesh.connectivityData.owner_[faceIndex]
					cellIndexN = mesh.connectivityData.neighbour_[faceIndex]
					
					# This holds only for a single component phi field!
					
					lower[cmpt][faceIndex] += min(-phi[0][faceIndex], 0)
					
					diag[cmpt][cellIndexP] += max(phi[0][faceIndex], 0)
					
					diag[cmpt][cellIndexN] += max(-phi[0][faceIndex], 0)
					
					upper[cmpt][faceIndex] += min(phi[0][faceIndex], 0)
			
			
			# Add boundary contributions. The loop goes through the boundary
			# patches and, depending on the boundaryType, defines the source and
			# diagonal contributions by calling the appropriate functions
			for boundaryPatch in mesh.connectivityData.boundary_:
				boundaryName = boundaryPatch[0]
				boundaryType = (psi.data.boundaryDefinition_[boundaryName])[0]
				
				boundaryContributionsFunctionString =						   \
					"self.boundaryContributions."			   \
					+ boundaryType											   \
					+ "(mesh, psi, boundaryPatch, source, diag, phi)"
				
				source, diag = eval(boundaryContributionsFunctionString)
					
			return source, lower, diag, upper


"""
// ************************************************************************* //
"""
