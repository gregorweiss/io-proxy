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
    Volume scalar field class.
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.volField.volField import *


class volScalarField(volField):
	
	fieldTypeName_ = "scalar"
	className_ = "volScalarField"
	noComponents_ = 1
	componentsNames_ = ("",)
	
	"""------------------------- Constructors -------------------------------"""
	
	# Constructor used to calculate the divergence of a volField
	@classmethod
	def div(self, field):
		
		mesh = field.data.mesh_
		
		fieldName = "div" + field.data.fieldName_
		
		boundaryDefinitionField = field.data.boundaryDefinition_
		boundaryDefinition = dict.fromkeys(boundaryDefinitionField)
		
		for fieldPatch in boundaryDefinitionField:
			 
			if (boundaryDefinitionField[fieldPatch][0] == 'empty'):
				boundaryDefinition[fieldPatch] = boundaryDefinitionField[fieldPatch]
			else:
				boundaryDefinition[fieldPatch] = 							   \
					('calculated', None)
		
		fieldClass = field.className_
		
		cellValuesFunc = eval("self.setCellValuesDiv" + fieldClass)
		boundaryValuesFunc = eval("self.setBoundaryValuesInitialDiv" + fieldClass)
		
		# Set cell values
		cellValues		= cellValuesFunc(field)
		
		# Set boundary values
		boundaryValues  = boundaryValuesFunc(cellValues, field)
		
		return self(mesh, fieldName, boundaryValues, cellValues, boundaryDefinition)
	
	
	"""------------------------ Initialization ------------------------------"""
	@classmethod
	def setCellValuesDivsurfaceScalarField(self, field):
		
		noComponents = 1
		
		mesh = field.data.mesh_
		
		nCells = mesh.geometryData.meshSize_
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		nFacesTot = np.size(mesh.connectivityData.owner_)
		
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		faceValues = field.data.faceValues_
		
		cellValues = np.zeros((noComponents, nCells), dtype = float)
		
		# Loop over internal faces
		for faceIndex in range(nInternalFaces):
			
			ownIndex = owner[faceIndex]
			neiIndex = neighbour[faceIndex]
			
			cellValues[0][ownIndex] += faceValues[0][faceIndex]
			cellValues[0][neiIndex] -= faceValues[0][faceIndex]
		
		# Loop over boundary faces
		for faceIndex in range(nInternalFaces, nFacesTot):
			
			boundaryCellIndex = owner[faceIndex]
			
			cellValues[0][boundaryCellIndex] += faceValues[0][faceIndex]
		
		cellValues[0] /= mesh.geometryData.V_
		
		return cellValues
		
	@classmethod
	def setBoundaryValuesInitialDivsurfaceScalarField(self, cellValues, field):
		
		noComponents = 1
		
		mesh = field.data.mesh_
		
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		nFacesTot = np.size(mesh.connectivityData.owner_)
		nFacesBound = nFacesTot - nInternalFaces
		
		owner = mesh.connectivityData.owner_
		
		boundaryValues = np.empty((noComponents, nFacesBound), dtype = float)
		
		for faceIndex in range(nInternalFaces, nFacesTot):
			
			boundaryFaceIndex = faceIndex - nInternalFaces
			
			boundaryCellIndex = owner[faceIndex]
			
			# Zero gradient extrapolation
			boundaryValues[0][boundaryFaceIndex] = cellValues[0][boundaryCellIndex]

		return boundaryValues
	
	
	"""----------------------- General functions ----------------------------"""
	
	
	# Write only one component
	def writeCellValues(self, file):
		
		nCells = self.data.mesh_.geometryData.meshSize_
		noComponents = self.noComponents_
		cellValues = self.data.cellValues_
		
		for cellIndex in range(nCells):

			# Write all components but the last with a trailing space symbol
			for cmpt in range(noComponents - 1):
				file.write(str(cellValues[cmpt][cellIndex]) + ' ')
			
			# Write the last component, close the bracket and insert new line
			file.write(str(cellValues[noComponents - 1][cellIndex]) + '\n')
			
	class writeBoundaryValuesFunctions:
		
		def fixedValue(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			
			fValue = (boundaryDefinition[boundaryName])[1]
			
			file.write(2*'\t' + 'value' + 2*'\t' + 'uniform ' +				   \
						   str(fValue[0]) + ';\n')
		
		def empty(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			None
			
		def fixedGradient(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			
			fGradient = (boundaryDefinition[boundaryName])[1]
			
			file.write(2*'\t' + 'gradient' + 2*'\t' + 'uniform ' +			   \
						   str(fGradient[0]) + ';\n')

"""
// ************************************************************************* //
"""
