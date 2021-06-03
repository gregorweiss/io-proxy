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
    Boundary contributions for the fvc gradient operator
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

class gradBoundaryContributions:
	
	class boundaryContributions:
		
		def fixedValue(mesh, psi, boundaryPatch, source, diag, field):
			
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			# Values of psi on the boundary
			fieldBoundaryValues = field.data.boundaryValues_
			
			C = mesh.geometryData.C_
			Cf = mesh.geometryData.Cf_
			Sf = mesh.geometryData.Sf_
			
			# Number of components of the field to make a gradient of
			noComponentsField = field.noComponents_
			
			# The tensor rank of the result is 1 higher than the field rank
			noComponentsGrad = field.noComponents_ * 3
			
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			owner = mesh.connectivityData.owner_
			
			# For all boundary faces
			for faceIndex in range(startFace, startFace + nFaces):
						
				boundaryCellIndex = owner[faceIndex]
						   
				boundaryIndex = faceIndex - nInternalFaces
				
				# For all grad components
				cmptGrad = 0
				# For all Sf components (3)
				for cmptSf in range(3):
					# For all psi components
					for cmptField in range(noComponentsField):
						# b_owner = Sf * field_b
						source[cmptGrad][boundaryCellIndex] +=				   \
							fieldBoundaryValues[cmptField][boundaryIndex] *    \
							Sf[faceIndex][cmptSf]
						
						# Increment gradient result component counter
						cmptGrad += 1
			
			return source, diag
			
			
		def empty(mesh, psi, boundaryPatch, source, diag, field):
			
			None
			
			return source, diag
			
		# Both the fixedGradient and fixedValue boundary contributions are
		# assembled using existing calculated boundary values
		#(psi.data.boundaryValues_), so it is essentially the same function:
		
		fixedGradient = fixedValue


"""
// ************************************************************************* //
"""
