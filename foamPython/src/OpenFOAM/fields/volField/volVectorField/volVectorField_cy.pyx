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
    Volume vector field class, with Cython.
    
\*---------------------------------------------------------------------------*/
"""

from src.OpenFOAM.fields.volField.volField import *
import numpy as np
cimport numpy as np

class volVectorField(volField):
	
	fieldTypeName_ = "vector"
	className_ = "volVectorField"
	noComponents_ = 3
	componentsNames_ = ("x", "y", "z")
	
	# Constructor used to calculate the gradient of a volField
	@classmethod
	def grad(self, field):
		
		mesh = field.data.mesh_
		
		fieldName = "grad" + field.data.fieldName_
		
		boundaryDefinitionField = field.data.boundaryDefinition_
		boundaryDefinition = dict.fromkeys(boundaryDefinitionField)
		
		for fieldPatch in boundaryDefinitionField:
			 
			if (boundaryDefinitionField[fieldPatch][0] == 'empty'):
				boundaryDefinition[fieldPatch] = boundaryDefinitionField[fieldPatch]
			else:
				boundaryDefinition[fieldPatch] = 							   \
					('calculated', None)
		
		# Set cell values
		cellValues		= self.setCellValuesGrad(field)
		
		# Set boundary values
		boundaryValues  = self.setBoundaryValuesInitialGrad(cellValues, field)
		
		return self(mesh, fieldName, boundaryValues, cellValues, boundaryDefinition)
	
	
	@classmethod
	def setCellValuesGrad(self, field):
		
		"""------------------- Cython declarations --------------------------"""
#		cdef int nInternalFaces
#		cdef int faceIndex
#		cdef int ownIndex
#		cdef int neiIndex
		
#		cdef np.ndarray[np.float_t, ndim=1] psi
#		cdef np.ndarray[np.float_t, ndim=1] lower
#		cdef np.ndarray[np.float_t, ndim=1] diag
#		cdef np.ndarray[np.float_t, ndim=1] upper
#		cdef np.ndarray[np.int_t, ndim=1] owner
#		cdef np.ndarray[np.int_t, ndim=1] neighbour
#		cdef np.ndarray[np.float_t, ndim=1] resultArray
		
		cdef int noComponentsField
		cdef int noComponentsGrad
		cdef int nCells
		cdef int nInternalFaces
		cdef int nFacesTot
		cdef int faceIndex
		cdef int ownIndex
		cdef int neiIndex
		cdef float absPf
		cdef float absNf
		cdef float f
		cdef int cmpt
		cdef int cmptGrad
		cdef int cmptSf
		cdef int cmptField
		cdef int boundaryCellIndex
		cdef int boundaryFaceIndex
		
		cdef np.ndarray[np.float_t, ndim=2] C
		cdef np.ndarray[np.float_t, ndim=2] Cf
		cdef np.ndarray[np.float_t, ndim=2] Sf
		cdef np.ndarray[np.float_t, ndim=1] V
		cdef np.ndarray[np.float_t, ndim=2] cellValues
		cdef np.ndarray[np.float_t, ndim=2] gradField
		cdef np.ndarray[np.int_t, ndim=1] owner
		cdef np.ndarray[np.int_t, ndim=1] neighbour
		cdef np.ndarray[np.float_t, ndim=1] interpCellValues
		cdef np.ndarray[np.float_t, ndim=2] boundaryValues
		
		"""------------------------------------------------------------------"""
		
		mesh = field.data.mesh_
		
		C = mesh.geometryData.C_
		Cf = mesh.geometryData.Cf_
		Sf = mesh.geometryData.Sf_
		V = mesh.geometryData.V_
		
		# Number of components of the field to make a gradient of
		noComponentsField = field.noComponents_
		
		# The tensor rank of the result is 1 higher than the field rank
		noComponentsGrad = field.noComponents_ * 3
		
		# Total number of cells
		nCells = mesh.read("geometryData.meshSize")
		# Number of internal faces
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		# Total number of faces
		nFacesTot = np.size(mesh.connectivityData.owner_)
		
		cellValues = field.data.cellValues_
		
		gradField = np.zeros((noComponentsGrad, nCells), dtype = float)
		
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
							
				interpCellValues[cmpt] = 									   \
					f * cellValues[cmpt][ownIndex] + 						   \
					(1 - f) * cellValues[cmpt][neiIndex]
			
			# For all grad components
			cmptGrad = 0
			# For all Sf components (3)
			for cmptSf in range(3):
				# For all field components
				for cmptField in range(noComponentsField):
					# int_V gradField dV = sum_f Sf * field_f > owner
					gradField[cmptGrad][ownIndex] +=						   \
						interpCellValues[cmptField] * 						   \
						Sf[faceIndex][cmptSf]
					
					# int_V gradField dV = - sum_f Sf * field_f > neighbour
					gradField[cmptGrad][neiIndex] -=						   \
						interpCellValues[cmptField] * 						   \
						Sf[faceIndex][cmptSf]
					
					# Increment gradient result component counter
					cmptGrad += 1
		
		boundaryValues = field.data.boundaryValues_
		
		# For all boundary faces
		for faceIndex in range(nInternalFaces, nFacesTot):
			
			boundaryCellIndex = owner[faceIndex]
			
			boundaryFaceIndex = faceIndex - nInternalFaces
			
			# For all grad components
			cmptGrad = 0
			# For all Sf components (3)
			for cmptSf in range(3):
				# For all field components
				for cmptField in range(noComponentsField):
					# int_V gradField dV = sum_f Sf * field_f
					gradField[cmptGrad][boundaryCellIndex] +=				   \
						boundaryValues[cmptField][boundaryFaceIndex] *		   \
						Sf[faceIndex][cmptSf]
						
					# Increment gradient result component counter
					cmptGrad += 1
			
		
		# Divide all cell values by the volume of the specific cell. The
		# previously calculated gradient values were actually volume integrals
		for cmptGrad in range(noComponentsGrad):
			gradField[cmptGrad] /= V
		
		
		return gradField
	
	
	@classmethod
	def setBoundaryValuesInitialGrad(self, cellValues, field):
		
		boundaryDefinition = field.data.boundaryDefinition_
		
		mesh = field.data.mesh_
		
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
				"self.setBoundaryValuesFunctionsGrad." + boundaryType)
			
			setBoundaryValuesFunc											   \
				(mesh, noComponents, cellValues, field, boundaryValues, boundaryPatch, boundaryDefinition)
		
		return boundaryValues
	
	
	class writeBoundaryValuesFunctions:
		
		def fixedValue(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			
			fValue = (boundaryDefinition[boundaryName])[1]
			
			file.write(2*'\t' + 'value' + 2*'\t' + 'uniform ')
			file.write('(')
			# X component
			file.write(str(fValue[0]) + ' ')
			# Y component
			file.write(str(fValue[1]) + ' ')
			# Z component
			file.write(str(fValue[2]) + ');\n')
		
		def empty(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			None
			
		def fixedGradient(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			
			fGradient = (boundaryDefinition[boundaryName])[1]
			
			file.write(2*'\t' + 'gradient' + 2*'\t' + 'uniform ')
			file.write('(')
			# X component
			file.write(str(fGradient[0]) + ' ')
			# Y component
			file.write(str(fGradient[1]) + ' ')
			# Z component
			file.write(str(fGradient[2]) + ');\n')
			
		def calculated(file, boundaryName, boundaryDefinition, mesh, boundaryValues):
			
			for i in mesh.connectivityData.boundary_:
				
				if (i[0] == boundaryName):
					boundaryPatch = i
					break
					
			nInternalFaces = np.size(mesh.connectivityData.neighbour_)
			
			nFaces = boundaryPatch[2]
			startFace = boundaryPatch[3]
			
			owner = mesh.connectivityData.owner_
			neighbour = mesh.connectivityData.neighbour_
			
			file.write(2*'\t' + 'value' + 2*'\t' + 'nonuniform List<vector>\n')
			file.write(str(nFaces) + '\n')
			file.write('(\n')
			
			for faceIndex in range(startFace, startFace + nFaces):
				
				boundaryFaceIndex = faceIndex - nInternalFaces
				
				file.write('(')
				# X component
				file.write(str(boundaryValues[0][boundaryFaceIndex]) + ' ')
				# Y component
				file.write(str(boundaryValues[1][boundaryFaceIndex]) + ' ')
				# Z component
				file.write(str(boundaryValues[2][boundaryFaceIndex]) + ')\n')
			
			file.write(');\n')
			

"""
// ************************************************************************* //
"""
