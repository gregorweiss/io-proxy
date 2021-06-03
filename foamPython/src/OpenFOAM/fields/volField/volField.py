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
    Generic volume field class.
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.fields import *
from src.OpenFOAM.fields.volField.volFieldBoundaryConditions import *


class volField(fields, volFieldBoundaryConditions):
	
	fieldTypeName_ = "fieldType"
	className_ = "volField"
	noComponents_ = 0
	componentsNames_ = None
	
	class volFieldData:
		
		def __init__(self):
			
			mesh_					= None # The mesh
			
			fieldName_				= None # The name of the field
			cellValues_ 			= None # Value of field in cell centres
			boundaryValues_			= None # Value of field on the boundary
										   # patches
			boundaryDefinition_		= None # User specified in the case setup
										   # file
	
	
	"""------------------------- Constructors -------------------------------"""	
	# Main constructor
	def __init__															   \
		(self, mesh, fieldName, boundaryValues, cellValues, boundaryDefinition):
		
		self.data = self.volFieldData()
		
		self.data.mesh_					= mesh
		self.data.fieldName_	 		= fieldName
		self.data.cellValues_	 		= cellValues
		self.data.boundaryValues_ 		= boundaryValues	
		self.data.boundaryDefinition_	= boundaryDefinition
		
		noComponents = self.noComponents_

	# Constructor which initializes the field to a specified value in all cells
	@classmethod
	def initialize(self, fieldName, mesh, boundaryDefinition, internalField):
		
		# Set cell values
		cellValues	 	= self.setCellValuesSpecified(mesh, internalField)
		
		# Set boundary values
		boundaryValues  = self.setBoundaryValuesInitial(cellValues, mesh, boundaryDefinition)
		
		return self(mesh, fieldName, boundaryValues, cellValues, boundaryDefinition)
	
	
	"""----------------------- Initialization functions ---------------------"""
	
	@classmethod
	def setCellValuesSpecified(self, mesh, internalField):
		
		noComponents = self.noComponents_
		
		nCells = mesh.read("geometryData.meshSize")
		
		# Each component (XYZ) is its own array
		cellValues = np.ones((noComponents, nCells), dtype = float)
		
		for cmpt in range(noComponents):
			cellValues[cmpt] *= internalField[cmpt]
		
		return cellValues
	
	@classmethod
	def setBoundaryValuesInitial(self, cellValues, mesh, boundaryDefinition):
		
		# Define the required parameters and set the boundary values
		nFacesBound = np.size(mesh.connectivityData.owner_)					   \
					  - np.size(mesh.connectivityData.neighbour_)
		
		noComponents = self.noComponents_
		
		# Initialize the boundary array
		boundaryValues = np.empty((noComponents, nFacesBound), dtype = float)
		
		# Loop over the boundaries, and assign corresponding values,
		# depending on the boundary definition array
		for boundaryPatch in mesh.connectivityData.boundary_:
			boundaryName = boundaryPatch[0]
			boundaryType = (boundaryDefinition[boundaryName])[0]
			
			setBoundaryValuesFunc = eval(									   \
				"self.setBoundaryValuesFunctions." + boundaryType)
			
			setBoundaryValuesFunc											   \
				(mesh, noComponents, cellValues, boundaryValues, boundaryPatch, boundaryDefinition)
		
		return boundaryValues


	"""------------------------- General functions -----------------------"""
	def setBoundaryValues(self, boundaryDefinition):
		
		mesh = self.data.mesh_
		
		boundaryValues = self.data.boundaryValues_
		
		cellValues = self.data.cellValues_
		
		noComponents = self.noComponents_
		
		# Loop over the boundaries, and assign corresponding values,
		# depending on the boundary definition array
		for boundaryPatch in mesh.connectivityData.boundary_:
			boundaryName = boundaryPatch[0]
			boundaryType = (boundaryDefinition[boundaryName])[0]
			
			setBoundaryValuesFunc = eval(									   \
				"self.setBoundaryValuesFunctions." + boundaryType)
			
			setBoundaryValuesFunc											   \
				(mesh, noComponents, cellValues, boundaryValues, boundaryPatch, boundaryDefinition)
	
	
	def relax(self, fieldOld, fvSolutionParameters):
		
		try:
			alpha = fvSolutionParameters['expUR']
		except:
			alpha = 1.0
		
		oldCellValues = fieldOld.data.cellValues_
		newCellValues = self.data.cellValues_
		
		oldBoundaryValues = fieldOld.data.boundaryValues_
		newBoundaryValues = self.data.boundaryValues_
		
		resultCellValues =													   \
			oldCellValues + alpha * (newCellValues - oldCellValues)
			
		resultBoundaryValues = 												   \
			oldBoundaryValues + alpha * (newBoundaryValues - oldBoundaryValues)
		
		self.data.cellValues_ = resultCellValues
		self.data.boundaryValues_ = resultBoundaryValues
	
	
	def setRefValue(self, fvSolution):
		
		boundaryDefinition = self.data.boundaryDefinition_
		
		# Loop over the boundaries and simply exit the function if there is a
		# fixedValue boundary
		for patch in boundaryDefinition:
			
			if ((boundaryDefinition[patch])[0] == 'fixedValue'):
				return
		
		noComponents = self.noComponents_
		
		# Read refCell and readValue, otherwise set them to zero
		try:
			refCell = fvSolution['refCell']
			refValue = fvSolution['refValue']
		except:
			refCell = 0
			refValue = np.zeros(noComponents, dtype = float)
		
		cellValues = self.data.cellValues_
		boundaryValues = self.data.boundaryValues_
		
		nCells = self.data.mesh_.geometryData.meshSize_
		
		for cmpt in range(noComponents):
			deltaValue = refValue[cmpt] - cellValues[cmpt][refCell]
			
			cellValues[cmpt] += deltaValue
			boundaryValues[cmpt] += deltaValue
			
	
	# Copy an existing field
	def copy(self):
		
		noComponents = self.noComponents_
		
		from src.OpenFOAM.fields.include import volScalarField, volVectorField
		
		resultClass = eval(self.className_)
		
		mesh = self.data.mesh_
		
		cellValues = np.copy(self.data.cellValues_)
		boundaryValues = np.copy(self.data.boundaryValues_)
		
		boundaryDefinition = self.data.boundaryDefinition_

		fieldName = self.data.fieldName_
		
		return resultClass(mesh, fieldName, boundaryValues, cellValues, boundaryDefinition)
	
	
	# Multiply all cell and boundary values with cell volume
	def integrateVolume(self):
		
		mesh = self.data.mesh_
		
		noComponents = self.noComponents_
		
		# Multiply cell values by volume
		for cmpt in range(noComponents):
			
			self.data.cellValues_[cmpt] *= mesh.geometryData.V_
		
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		nFacesTot = np.size(mesh.connectivityData.owner_)
		
		owner = mesh.connectivityData.owner_
		
		# Multiply boundary values by volume of the owner
		for faceIndex in range(nInternalFaces, nFacesTot):
			
			boundaryCellIndex = owner[faceIndex]
			
			boundaryFaceIndex = faceIndex - nInternalFaces
			
			for cmpt in range(noComponents):
				
				self.data.boundaryValues_[cmpt][boundaryFaceIndex] *=		   \
					mesh.geometryData.V_[boundaryCellIndex]
	
	
	# Divide all cell and boundary values by cell volume
	def divideVolume(self):
		
		mesh = self.data.mesh_
		
		noComponents = self.noComponents_
		
		# Divide cell values by volume
		for cmpt in range(noComponents):
			
			self.data.cellValues_[cmpt] /= mesh.geometryData.V_
		
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		nFacesTot = np.size(mesh.connectivityData.owner_)
		
		owner = mesh.connectivityData.owner_
		
		# Divide boundary values by volume of the owner
		for faceIndex in range(nInternalFaces, nFacesTot):
			
			boundaryCellIndex = owner[faceIndex]
			
			boundaryFaceIndex = faceIndex - nInternalFaces
			
			for cmpt in range(noComponents):
				
				self.data.boundaryValues_[cmpt][boundaryFaceIndex] /=		   \
					mesh.geometryData.V_[boundaryCellIndex]
	
	
	# Write field values into the corresponding time folder
	# If a file with the same name exists in the time folder, it is overwritten
	def write(self, writeTime):
		
		# Local variable names
		fieldTypeName = self.fieldTypeName_
		className = self.className_
		noComponents = self.noComponents_
		
		fieldName = self.data.fieldName_
		
		nCells = self.data.mesh_.geometryData.meshSize_
		
		# Create time directory if it does not exist
		if os.path.exists(str(writeTime)):
			None
		else:
			os.mkdir(str(writeTime))
		
		# Open the new file
		file = open(str(writeTime) +'/' + fieldName, 'w')
		
		# Write the header and beginning of file
		file.write(HEADER)
		
		file.write('    class       ' + className +';\n')
		file.write('    location    ' + '"' + str(writeTime) + '"' + ';\n')
		file.write('    object      ' + fieldName + ';\n}\n')
		file.write('// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n')
		file.write('dimensions      [0 0 0 0 0 0 0];\n\n')
		file.write('internalField   nonuniform List<' + fieldTypeName + '>\n')
		
		file.write(str(nCells) + '\n(\n')
		
		# This function is in the individual class (volScalarField, etc.)
		self.writeCellValues(file)
			
		file.write(');\n\n')
		file.write('boundaryField\n{\n')
		
		self.writeBoundaryValues(file)
			
		file.write('}\n\n')
		file.write('// ************************************************************************* //')
	
	# Write all components in the appropriate format
	def writeCellValues(self, file):
		
		nCells = self.data.mesh_.geometryData.meshSize_
		noComponents = self.noComponents_
		cellValues = self.data.cellValues_
		
		for cellIndex in range(nCells):

			file.write('(')
			
			# Write all components but the last with a trailing space symbol
			for cmpt in range(noComponents - 1):
				file.write(str(cellValues[cmpt][cellIndex]) + ' ')
			
			# Write the last component, close the bracket and insert new line
			file.write(str(cellValues[noComponents - 1][cellIndex]) + ')\n')
	
	# Used to call the appropriate function in the individual class
	# (volScalarField, etc.)
	def writeBoundaryValues(self, file):
		
		boundaryDefinition = self.data.boundaryDefinition_
		
		boundaryValues = self.data.boundaryValues_
		
		mesh = self.data.mesh_
		
		# Write boundary values loop
		for boundaryName in boundaryDefinition:
			
			boundaryType = (boundaryDefinition[boundaryName])[0]
			
			file.write('\t' + boundaryName + '\n')
			file.write('\t' + '{\n')
			file.write(2*'\t' + 'type' + 2*'\t' +							   \
						(boundaryDefinition[boundaryName])[0] + ';\n')
			
			# Assign and call the appropriate function, depending on the
			# boundaryType
			writeBoundaryValuesFunc = eval(									   \
				"self.writeBoundaryValuesFunctions." + 					   \
				boundaryType)
				
			writeBoundaryValuesFunc(file, boundaryName, boundaryDefinition, mesh, boundaryValues)

			file.write('\t' + '}\n\n')
	
	
	"""------------------------ Defining operators --------------------------"""
	
	def __sub__(self, other): # - operator
					
		noComponents = self.noComponents_
		
		if (self.noComponents_ != other.noComponents_):
			raise RuntimeError('Only volume fields with the same number of' +  \
				' components can be subtracted!')
				
		from src.OpenFOAM.fields.include import volScalarField, volVectorField
		
		if (noComponents == 1):
			resultClass = volScalarField
		elif (noComponents == 3):
			resultClass = volVectorField
		
		mesh = self.data.mesh_
		
		cellValues = self.data.cellValues_ - other.data.cellValues_
		boundaryValues = self.data.boundaryValues_ - other.data.boundaryValues_
		
		selfBoundaryDefinition = self.data.boundaryDefinition_
		resultBoundaryDefinition = dict.fromkeys(selfBoundaryDefinition)
		
		for selfPatch in selfBoundaryDefinition:
			 
			if (selfBoundaryDefinition[selfPatch][0] == 'empty'):
				resultBoundaryDefinition[selfPatch] = selfBoundaryDefinition[selfPatch]
			else:
				resultBoundaryDefinition[selfPatch] = 							   \
					('calculated', None)
		
		selfName = self.data.fieldName_
		otherName = other.data.fieldName_
		fieldName = selfName + "-" + otherName
		
		return resultClass(mesh, fieldName, boundaryValues, cellValues, resultBoundaryDefinition)
		
		
	def __mul__(self, other): # * operator
		
		from src.OpenFOAM.fields.include import volScalarField, volVectorField
		
		noComponentsSelf = self.noComponents_
		noComponentsOther = other.noComponents_
		
		# The number of components of the resulting field
		noComponents = max(noComponentsSelf, noComponentsOther)
		
		# Choose the appropriate resulting class, depending on the class of the field
		if (noComponents == 1):
			resultClass = volScalarField
		elif (noComponents == 3):
			resultClass = volVectorField
		else:
			raise RuntimeError("Error in multiplication of vector fields!")
		
		selfName = self.data.fieldName_
		otherName = other.data.fieldName_
		resultName = selfName + "*" + otherName
		
		mesh = self.data.mesh_

		resultInternalField = np.zeros(noComponents, dtype=float)
		
		selfBoundaryDefinition = self.data.boundaryDefinition_
		resultBoundaryDefinition = dict.fromkeys(selfBoundaryDefinition)
		
		for selfPatch in selfBoundaryDefinition:
			 
			if (selfBoundaryDefinition[selfPatch][0] == 'empty'):
				resultBoundaryDefinition[selfPatch] = selfBoundaryDefinition[selfPatch]
			else:
				resultBoundaryDefinition[selfPatch] = 							   \
					('calculated', None)
		
		resultVolField = resultClass.initialize(							   \
			resultName, mesh, resultBoundaryDefinition, resultInternalField)
		
		nCells = self.data.mesh_.geometryData.meshSize_
		
		nBoundaryFaces = 													   \
			np.size(mesh.connectivityData.owner_) - np.size(mesh.connectivityData.neighbour_)
		
		resultCellValues = np.empty((noComponents, nCells), dtype = float)
		resultBoundaryValues = np.empty((noComponents, nBoundaryFaces), dtype = float)
		
		selfCellValues = self.data.cellValues_
		selfBoundaryValues = self.data.boundaryValues_
		otherCellValues = other.data.cellValues_
		otherBoundaryValues = other.data.boundaryValues_
		
		resultCmpt = 0
		for selfCmpt in range(noComponentsSelf):
			for otherCmpt in range(noComponentsOther):
				
				resultCellValues[resultCmpt] = selfCellValues[selfCmpt] * otherCellValues[otherCmpt]
				resultBoundaryValues[resultCmpt] = selfBoundaryValues[selfCmpt] * otherBoundaryValues[otherCmpt]
				
				resultCmpt += 1
				
		resultVolField.data.cellValues_ = resultCellValues
		resultVolField.data.boundaryValues_ = resultBoundaryValues
		
		return resultVolField
		

"""
// ************************************************************************* //
"""
