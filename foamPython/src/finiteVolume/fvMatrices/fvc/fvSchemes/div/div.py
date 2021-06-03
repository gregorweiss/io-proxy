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
    Class which contains the divergence contributions for fvc
    
\*---------------------------------------------------------------------------*/
"""
import numpy as np

class divClass:
	
	class div():
		
		@classmethod
		def linear(self, mesh, psi, fvVariables):
			
			field = fvVariables[0]
			
			C = mesh.geometryData.C_
			Cf = mesh.geometryData.Cf_
			Sf = mesh.geometryData.Sf_
			
			# Number of components of the field to make a gradient of
			noComponentsField = field.noComponents_
			
			# The tensor rank of the result is 1 lower than the field rank
			noComponentsDiv = int(field.noComponents_ / 3 + 1e-10)
			
			if (noComponentsDiv != 1):
				raise RuntimeError("Wrong or not yet implemented field type for div!")
			
			# Total number of cells
			nCells = mesh.read("geometryData.meshSize")
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			nFacesTot = np.size(mesh.connectivityData.owner_)
			
			cellValues = field.data.cellValues_
			boundaryValues = field.data.boundaryValues_
			
			source = np.zeros((noComponentsDiv, nCells), dtype = float)
			
			lower = np.zeros((noComponentsDiv, nInternalFaces), dtype = float)
			diag = np.zeros((noComponentsDiv, nCells), dtype = float)
			upper = np.zeros((noComponentsDiv, nInternalFaces), dtype = float)
			
			owner = mesh.connectivityData.owner_
			neighbour = mesh.connectivityData.neighbour_
			
			faceValue = np.empty(noComponentsField, dtype = float)
			
			# For all internal faces
			for faceIndex in range(nInternalFaces):
					
				ownIndex = owner[faceIndex]
				neiIndex = neighbour[faceIndex]
				
				absPf = np.linalg.norm(Cf[faceIndex] - C[ownIndex])
				absNf = np.linalg.norm(Cf[faceIndex] - C[neiIndex])
				
				f = absNf / (absPf + absNf)
				
				# For all field components, interpolate value to face
				for cmptField in range(noComponentsField):
								
					faceValue[cmptField] = 									   \
						f * cellValues[cmptField][ownIndex] + 				   \
						(1 - f) * cellValues[cmptField][neiIndex]
				
				# Works only for a scalar equation (div of a vector field)!!!
				
				# Div component
				cmptDiv = 0
				
				source[cmptDiv][ownIndex] += 								   \
					np.dot(faceValue, Sf[faceIndex])
					
				source[cmptDiv][neiIndex] -= 								   \
					np.dot(faceValue, Sf[faceIndex])
					
			for faceIndex in range(nInternalFaces, nFacesTot):
				
				boundaryFaceIndex = faceIndex - nInternalFaces
				
				boundaryCellIndex = owner[faceIndex]
				
				# Works only for a scalar equation (div of a vector field)!!!
				
				# Div component
				cmptDiv = 0
				
				for cmptField in range(noComponentsField):
					faceValue[cmptField] = boundaryValues[cmptField][boundaryFaceIndex]
				
				source[cmptDiv][boundaryCellIndex] +=						   \
					np.dot(faceValue, Sf[faceIndex])
					
			return source, lower, diag, upper


"""
// ************************************************************************* //
"""
