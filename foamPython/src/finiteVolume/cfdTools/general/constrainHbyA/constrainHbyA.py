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
    Function to constrain HbyA in pressure-velocity coupling
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

def constrainHbyA(HbyA, U):
	
	boundaryDefinition_U = U.data.boundaryDefinition_
	
	# """---------------- Go through the faces "by hand" ----------------------"""
	
	# meshBoundary = U.data.mesh_.connectivityData.boundary_
	
	# noComponents = HbyA.noComponents_
	
	# nInternalFaces = np.size(U.data.mesh_.connectivityData.neighbour_)
	
	# for patch in meshBoundary:
		
		# boundaryName = patch[0]
		
		# if ( (boundaryDefinition_p[boundaryName])[0] == 'fixedGradient' ):
			
			# startFace = patch[3]
			# nFaces = patch[2]
			
			# for cmpt in range(noComponents):
				
				# for faceIndex in range(startFace, startFace + nFaces):
					
					# boundaryFaceIndex = faceIndex - nInternalFaces
					
					# HbyA.data.boundaryValues_[cmpt][boundaryFaceIndex] =	   \
						# U.data.boundaryValues_[cmpt][boundaryFaceIndex]

	# """----------------------------------------------------------------------"""
	
	"""---------------- Set velocity boundary conditions --------------------"""
	
	HbyA.updateBoundaryDefinition(boundaryDefinition_U)
	HbyA.setBoundaryValues(boundaryDefinition_U)
