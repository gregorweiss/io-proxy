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
    Surface scalar field class
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.surfaceField.surfaceField import *


class surfaceScalarField(surfaceField):
	
	fieldTypeName_ = "scalar"
	className_ = "surfaceScalarField"
	noComponents_ = 1
	componentsNames_ = ("",)
	
	# Constructor which initializes the flux surface field using a specified
	# volVectorField
	@classmethod
	def flux(self, field):
		
		mesh = field.data.mesh_
		
		fieldName = "phi" + field.data.fieldName_
		
		boundaryDefinitionField = field.data.boundaryDefinition_
		boundaryDefinition = dict.fromkeys(boundaryDefinitionField)
		
		for patch in boundaryDefinitionField:
			 
			if (boundaryDefinitionField[patch][0] == 'empty'):
				boundaryDefinition[patch] = boundaryDefinitionField[patch]
			else:
				boundaryDefinition[patch] = 								   \
					('calculated', None)
		
		faceValues	 	= self.interpolateFlux(field)
		
		return self(mesh, fieldName, faceValues, boundaryDefinition)
		
		
	def updateFlux(self, volVectorField):
				
		self.data.faceValues_ = self.interpolateFlux(volVectorField)
	
	
	@classmethod
	def interpolateFlux(self, volVectorField):
		
		mesh = volVectorField.data.mesh_
		
		noComponents = volVectorField.noComponents_
		
		# Total number of faces
		nFacesTot = np.size(mesh.connectivityData.owner_)
		# Number of internal faces
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		C = mesh.geometryData.C_
		Cf = mesh.geometryData.Cf_
		Sf = mesh.geometryData.Sf_
		
		cellValues = volVectorField.data.cellValues_
		boundaryDefinition = volVectorField.data.boundaryDefinition_
		
		faceValues = np.empty((1, nFacesTot), dtype = float)
		
		interpVolVectorField = np.empty(noComponents, dtype = float)
		
		# For all internal faces
		for faceIndex in range(nInternalFaces):
				
			ownIndex = owner[faceIndex]
			neiIndex = neighbour[faceIndex]
			
			absPf = np.linalg.norm(Cf[faceIndex] - C[ownIndex])
			absNf = np.linalg.norm(Cf[faceIndex] - C[neiIndex])
			
			f = absNf / (absPf + absNf)
			
			# For all vector components, interpolate to face
			for cmpt in range(noComponents):			
				interpVolVectorField[cmpt] = 								   \
					f * cellValues[cmpt][ownIndex] + 						   \
					(1 - f) * cellValues[cmpt][neiIndex]
					
			# Result is the dot product of face area vector and the interpolated
			# volVectorField
			faceValues[0][faceIndex] =										   \
				np.dot(Sf[faceIndex], interpVolVectorField)
			
		
		# For all boundary faces
		
		volVectorFieldBoundaryValues = volVectorField.data.boundaryValues_
		
		for faceIndex in range(nInternalFaces, nFacesTot):
			
			boundaryIndex = faceIndex - nInternalFaces
			
			for cmpt in range(noComponents):
				interpVolVectorField[cmpt] = volVectorFieldBoundaryValues[cmpt][boundaryIndex]
			
			faceValues[0][faceIndex] =										   \
				np.dot(Sf[faceIndex], interpVolVectorField)
			
		return faceValues

"""
// ************************************************************************* //
"""
