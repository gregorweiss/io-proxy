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
    volField class boundary conditions functions
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

class volFieldBoundaryConditions:
	
	# Functions used to set the appropriate values in the boundaryValues array.
	# This is later used for writing results, calculating gradient, etc.
	class setBoundaryValuesFunctions:
		
		def fixedValue														   \
		   (mesh, noComponents, cellValues, boundaryValues, boundaryPatch, boundaryDefinition):

			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			# Boundary properties
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			# Start and end indices inside the boundaryValues array
			boundaryValuesStart = startFace - nInternalFaces
			boundaryValuesEnd = boundaryValuesStart + nFaces
			
			owner = mesh.connectivityData.owner_
			
			specifiedFixedValue = (boundaryDefinition[boundaryName])[1]
			
			# Assign the appropriate values to the appropriate part of the 
			# boundaryValues array
			for cmpt in range(noComponents):
				
				for faceIndex in range(startFace, startFace + nFaces):
					
					boundaryFaceIndex = faceIndex - nInternalFaces
					
					boundaryCellIndex = owner[faceIndex]
					
					boundaryValues[cmpt][boundaryFaceIndex]					   \
						= specifiedFixedValue[cmpt]
					
					
		def empty															   \
		   (mesh, noComponents, cellValues, boundaryValues, boundaryPatch, boundaryDefinition):
			   
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			owner = mesh.connectivityData.owner_
			
			# Boundary properties
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			for cmpt in range(noComponents):
				for faceIndex in range(startFace, startFace + nFaces):
					
					boundaryFaceIndex = faceIndex - nInternalFaces
					
					boundaryCellIndex = owner[faceIndex]
					
					boundaryValues[cmpt][boundaryFaceIndex]					   \
						= cellValues[cmpt][boundaryCellIndex]
						
		def fixedGradient													   \
		   (mesh, noComponents, cellValues, boundaryValues, boundaryPatch, boundaryDefinition):
			   
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			owner = mesh.connectivityData.owner_
			
			# Boundary properties
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			# Cell centres
			C = mesh.geometryData.C_
			# Face centres
			Cf = mesh.geometryData.Cf_
			
			specifiedFixedGradient = (boundaryDefinition[boundaryName])[1]
			
			# Assign the appropriate values to the appropriate part of the 
			# boundaryValues array
			for cmpt in range(noComponents):
				for faceIndex in range(startFace, startFace + nFaces):
					
					boundaryFaceIndex = faceIndex - nInternalFaces
					
					boundaryCellIndex = owner[faceIndex]
					
					boundaryValues[cmpt][boundaryFaceIndex]					   \
						= cellValues[cmpt][boundaryCellIndex] + 			   \
						specifiedFixedGradient[cmpt] *						   \
						np.linalg.norm(Cf[faceIndex] - C[boundaryCellIndex])
						
		def calculated														   \
			(mesh, noComponents, cellValues, boundaryValues, boundaryPatch, boundaryDefinition):
			
			None
	
	# Functions used to set the appropriate values in the boundaryValues array
	# when calculating the gradient of a field.
	class setBoundaryValuesFunctionsGrad:
		
		def fixedValue														   \
		   (mesh, noComponentsGrad, cellValues, field, boundaryValues, boundaryPatch, boundaryDefinition):
			   
			C = mesh.geometryData.C_
			Cf = mesh.geometryData.Cf_
			
			noComponentsField = field.noComponents_
			
			fieldCellValues = field.data.cellValues_
			fieldBoundaryValues = field.data.boundaryValues_
			
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			# Boundary properties
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			owner = mesh.connectivityData.owner_
			
			# For all boundary faces, directly evaluate the gradient using given 
			# field values
			for faceIndex in range(startFace, startFace + nFaces):
				
				boundaryCellIndex = owner[faceIndex]
				
				boundaryFaceIndex = faceIndex - nInternalFaces
				
				# For all grad components
				cmptGrad = 0
				# For all field components
				for cmptField in range(noComponentsField):
					# For all distance vector components (3)
					for cmptD in range(3):
						
						cellCentre = C[boundaryCellIndex]
						faceCentre = Cf[boundaryCellIndex]
						
						gradFaceCell =										   \
							(fieldCellValues[cmptField][boundaryCellIndex] - fieldBoundaryValues[cmptField][boundaryFaceIndex])
							
						distFaceCell = 										   \
							(C[boundaryCellIndex][cmptD] - Cf[boundaryFaceIndex][cmptD])
						
						# If the face centre and the cell centre are at the same
						# location in a given coordinate axis, copy the value of
						# the cell gradient to the face
						if (abs(distFaceCell) < 1e-10):
							boundaryValues[cmptGrad][boundaryFaceIndex] =	   \
								cellValues[cmptGrad][boundaryCellIndex]
						else:
							boundaryValues[cmptGrad][boundaryFaceIndex] =	   \
							 (fieldCellValues[cmptField][boundaryCellIndex] - fieldBoundaryValues[cmptField][boundaryFaceIndex]) / \
							 (C[boundaryCellIndex][cmptD] - Cf[boundaryFaceIndex][cmptD])
						
						cmptGrad += 1
		
		
		# All face gradients are calculated using the pre-existing field values
		# in cells and boundary faces. For that reason, boundary calculation
		# functions are the same for all boundary conditions:
		
		empty = fixedValue
		
		fixedGradient = fixedValue

"""
// ************************************************************************* //
"""
