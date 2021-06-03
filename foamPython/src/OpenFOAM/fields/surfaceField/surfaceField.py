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
    Generic surface field class.
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.fields import *
from src.OpenFOAM.fields.surfaceField.surfaceFieldBoundaryConditions import *


class surfaceField(fields, surfaceFieldBoundaryConditions):
	
	fieldTypeName_ = "fieldType"
	className_ = "surfaceField"
	noComponents_ = 0
	componentsNames_ = None
	
	class surfaceFieldData:
		
		def __init__(self):
			
			mesh_					= None # The mesh
			
			fieldName_				= None # The name of the field
			
			faceValues_ 			= None # Value of field in face centres

			boundaryDefinition_		= None # User specified in the case setup
										   # file

	
	"""------------------------- Constructors -------------------------------"""	
	# Main constructor
	def __init__															   \
		(self, mesh, fieldName, faceValues, boundaryDefinition):
		
		self.data = self.surfaceFieldData()
		
		self.data.mesh_					= mesh
		
		self.data.fieldName_	 		= fieldName
		
		self.data.faceValues_	 		= faceValues

		self.data.boundaryDefinition_	= boundaryDefinition
	
	# Good only for surfaceScalarField
	def write(self, writeTime):
		
		# Local variable names
		fieldTypeName = self.fieldTypeName_
		className = self.className_
		noComponents = self.noComponents_
		
		fieldName = self.data.fieldName_
		
		nInternalFaces = np.size(self.data.mesh_.connectivityData.neighbour_)
		
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
		
		file.write(str(nInternalFaces) + '\n(\n')
		
		# This function is in the individual class (volScalarField, etc.)
		self.writeInternalFaceValues(file)
			
		file.write(');\n\n')
		file.write('boundaryField\n{\n')
		
		self.writeBoundaryValues(file)
			
		file.write('}\n\n')
		file.write('// ************************************************************************* //')
	
	# Good only for surfaceScalarField
	def writeInternalFaceValues(self, file):
		
		nInternalFaces = np.size(self.data.mesh_.connectivityData.neighbour_)
		noComponents = self.noComponents_
		faceValues = self.data.faceValues_
		
		for faceIndex in range(nInternalFaces):

			file.write(str(faceValues[0][faceIndex]) + '\n')
	
	# Good only for surfaceScalarField
	def writeBoundaryValues(self, file):
		
		boundaryDefinition = self.data.boundaryDefinition_
		
		faceValues = self.data.faceValues_
		
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
				
			writeBoundaryValuesFunc(file, boundaryName, boundaryDefinition, mesh, faceValues)

			file.write('\t' + '}\n\n')
	
	class writeBoundaryValuesFunctions:
		
		def calculated(file, boundaryName, boundaryDefinition, mesh, faceValues):
			
			for i in mesh.connectivityData.boundary_:
				
				if (i[0] == boundaryName):
					boundaryPatch = i
					break
			
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			file.write(2*'\t' + 'value' + 2*'\t' + 'nonuniform List<scalar>\n')
			file.write(str(nFaces) + '\n')
			file.write('(\n')
			
			for faceIndex in range(startFace, startFace + nFaces):
				
				file.write(str(faceValues[0][faceIndex]) + '\n')
				
			file.write(');\n')
		
		def empty(file, boundaryName, boundaryDefinition, mesh, faceValues):
			
			None

	
	"""------------------------ Defining operators --------------------------"""
	
	def __sub__(self, other): # - operator
		
		noComponents = self.noComponents_
		
		if (self.noComponents_ != self.noComponents_):
			raise RuntimeError('Only surface fields with the same number of' + \
				' components can be subtracted!')
				
		from src.OpenFOAM.fields.include import surfaceScalarField, surfaceVectorField

		if (noComponents == 1):
			resultClass = surfaceScalarField
		elif (noComponents == 3):
			resultClass = surfaceVectorField
		
		mesh = self.data.mesh_
		
		resultFaceValues = self.data.faceValues_ - other.data.faceValues_
		
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
		resultName = selfName + "-" + otherName
		
		return resultClass(mesh, resultName, resultFaceValues, resultBoundaryDefinition)

"""
// ************************************************************************* //
"""
