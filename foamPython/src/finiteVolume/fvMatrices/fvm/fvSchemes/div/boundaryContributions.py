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
    Boundary contributions for the fvm divergence operator
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

class divBoundaryContributions:
	
	class boundaryContributions:
		
		def fixedValue(mesh, psi, boundaryPatch, source, diag, phi):
			
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			# Values of psi on the boundary
			psiBoundaryValues = psi.data.boundaryValues_
			
			noComponents = psi.noComponents_
			
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			for cmpt in range(noComponents):
			
				for faceIndex in range(startFace, startFace + nFaces):
					
					boundaryCellIndex = mesh.connectivityData.owner_[faceIndex]
						   
					boundaryIndex = faceIndex - nInternalFaces
					
					# This holds only for a single component phi field!
					
					# b = F * psi_b
					source[cmpt][boundaryCellIndex] += 				  		   \
						phi[0][faceIndex] * psiBoundaryValues[cmpt][boundaryIndex]
			
			return source, diag
			
			
		def empty(mesh, psi, boundaryPatch, source, diag, phi):
			
			None
			
			return source, diag
			
			
		def fixedGradient(mesh, psi, boundaryPatch, source, diag, phi):
			
			boundaryName = boundaryPatch[0]
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			# Cell centres
			C = mesh.geometryData.C_
			# Face centres
			Cf = mesh.geometryData.Cf_
			
			# Number of internal faces
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			# The value of the gradient
			gradPsi_b = (psi.data.boundaryDefinition_[boundaryName])[1]
			
			noComponents = psi.noComponents_
			
			for cmpt in range(noComponents):
			
				for faceIndex in range(startFace, startFace + nFaces):
					
					boundaryCellIndex = mesh.connectivityData.owner_[faceIndex]
						   
					boundaryIndex = faceIndex - nInternalFaces
					
					# This holds only for a single component phi field!
					
					# ap = F
					diag[cmpt][boundaryCellIndex] += 						   \
						 phi[0][faceIndex]
					
					# b = F * |d_b| * gradPsi_b
					
					d_b = np.linalg.norm(Cf[faceIndex] - C[boundaryCellIndex])
					
					source[cmpt][boundaryCellIndex] += 				 	 	   \
						 phi[0][faceIndex] * d_b * gradPsi_b[cmpt]
						
			
			return source, diag

"""
// ************************************************************************* //
"""
