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
    Class which contains the gradient contributions for fvc
    
\*---------------------------------------------------------------------------*/
"""

from src.finiteVolume.fvMatrices.fvc.fvSchemes.grad.boundaryContributions import *

class gradClass:
	
	class grad(gradBoundaryContributions):
		
		@classmethod
		def linear(self, mesh, psi, fvVariables):
			
			field = fvVariables[0]
			
			C = mesh.geometryData.C_
			Cf = mesh.geometryData.Cf_
			Sf = mesh.geometryData.Sf_
			
			# Number of components of the field to make a gradient of
			noComponentsField = field.noComponents_
			
			# The tensor rank of the result is 1 higher than the field rank
			noComponentsGrad = field.noComponents_ * 3
			
			# Total number of cells
			nCells = mesh.read("geometryData.meshSize")
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			cellValues = field.data.cellValues_
			
			source = np.zeros((noComponentsGrad, nCells), dtype = float)
			
			lower = np.zeros((noComponentsGrad, nInternalFaces), dtype = float)
			diag = np.zeros((noComponentsGrad, nCells), dtype = float)
			upper = np.zeros((noComponentsGrad, nInternalFaces), dtype = float)
			
			owner = mesh.connectivityData.owner_
			neighbour = mesh.connectivityData.neighbour_
			
			interpCellValues = np.empty(noComponentsField, dtype = float)
			
			# For all internal faces
			for faceIndex in range(nInternalFaces):
					
				ownIndex = owner[faceIndex]
				neiIndex = neighbour[faceIndex]
				
				absPf = np.linalg.norm(Cf[faceIndex] - C[ownIndex])
				absNf = np.linalg.norm(Cf[faceIndex] - C[neiIndex])
				
				f = absNf / (absPf + absNf)
				
				# For all field components, interpolate value to face
				for cmpt in range(noComponentsField):
								
					interpCellValues[cmpt] = 								   \
						f * cellValues[cmpt][ownIndex] + 					   \
						(1 - f) * cellValues[cmpt][neiIndex]
				
				# For all grad components
				cmptGrad = 0
				# For all Sf components (3)
				for cmptSf in range(3):
					# For all field components
					for cmptField in range(noComponentsField):
						# b_owner = Sf * field_f
						source[cmptGrad][ownIndex] +=						   \
							interpCellValues[cmptField] * 					   \
							Sf[faceIndex][cmptSf]
						
						# b_neighbour = - Sf * field_f
						source[cmptGrad][neiIndex] -=						   \
							interpCellValues[cmptField] * 					   \
							Sf[faceIndex][cmptSf]
						
						# Increment gradient result component counter
						cmptGrad += 1
			
			# Add boundary contributions. The loop goes through the boundary
			# patches and, depending on the boundaryType, defines the source and
			# diagonal contributions by calling the appropriate functions
			for boundaryPatch in mesh.connectivityData.boundary_:
				boundaryName = boundaryPatch[0]
				boundaryType = (field.data.boundaryDefinition_[boundaryName])[0]
				
				boundaryContributionsFunctionString =						   \
					"self.boundaryContributions."			   \
					+ boundaryType											   \
					+ "(mesh, psi, boundaryPatch, source, diag, field)"
				
				source, diag = eval(boundaryContributionsFunctionString)
			
			return source, lower, diag, upper


"""
// ************************************************************************* //
"""
