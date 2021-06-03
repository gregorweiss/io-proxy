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
    A special matrix type and solver, designed for finite volume
    solutions of scalar equations.
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np

from src.finiteVolume.fvSolution.include import *
from src.OpenFOAM.fields.include import *

class fvMatrix():
	
	solversDict_ =															   \
	{																		   \
		"GaussSeidel" : GaussSeidel.solveMatrix,							   \
		"PointJacobi" : PointJacobi.solveMatrix								   \
	}
	
	class fvMatrixData:
		def __init__(self):
			
			mesh_		= None # fvMesh object
			
			psi_		= None # Field to be solved
			
			source_		= None # Source term
			
			# Matrix coefficients
			lower_		= None
			diag_		= None
			upper_		= None
		
			rowStart_	= None # Used to enable looping over cells/matrix rows

	"""------------------------- Constructors -------------------------------"""	
	# Main constructor
	def __init__(self, psi, source, lower, diag, upper, rowStart):
		
		self.data = self.fvMatrixData()
		
		self.data.mesh_ 				= psi.data.mesh_
		
		self.data.psi_					= psi
		
		self.data.source_				= source
		
		self.data.lower_				= lower
		self.data.diag_					= diag
		self.data.upper_				= upper
		self.data.rowStart_				= rowStart
	
	# Constructor to define the fvMatrix using string evaluation and a 
	# dictionary of operators which have been implemented within fvm or fvc
	@classmethod
	def construct(self, psi, operator, scheme, fvVariables):
		
		rowStart						= self.defineRowStart(psi.data.mesh_)
		
		source, lower, diag, upper		=									   \
			self.defineMatrix(psi, operator, scheme, fvVariables)
			
		return self(psi, source, lower, diag, upper, rowStart)


	"""------------- Matrix components definition functions -----------------"""
	
	# Calls the right solver based on the "solver" input specified for the field
	
	def solve(self, fvSolutionParameters):
		
		solversDict = self.solversDict_
		
		try:
			solver = fvSolutionParameters['solver']
		except:
			raise RuntimeError(												   \
				"No solver specified for field " + equationName)
		
		# Assign the appropriate solver to the solverFunction variable
		try:
			solverFunction = solversDict[solver]
		except:
			raise RuntimeError(												   \
				"Solver " + solver + " has not been implemented!")
		
		print() # Output new line
		
		# Run the appropriate solver
		solverFunction(self, fvSolutionParameters)
	
	
	def relax(self, fvSolutionParameters):
		
		alpha = 1.0
		
		try:
			alpha = fvSolutionParameters['impUR']
		except:
			None
		
		lower = self.data.lower_
		diag = self.data.diag_
		upper = self.data.upper_
		source = self.data.source_
		
		# Store the current unrelaxed diagonal for use in updating the source
		diagOld = np.copy(diag)
		
		mesh = self.data.psi_.data.mesh_
		nCells = mesh.geometryData.meshSize_
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		noComponents = self.data.psi_.noComponents_
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		# # sumMagOffDiag
		# sumOff = np.zeros((noComponents, nCells), dtype=float)

		# for cmpt in range(noComponents):
			
			# for faceIndex in range(nInternalFaces):
				
				# ownIndex = owner[faceIndex]
				# neiIndex = neighbour[faceIndex]
				
				# sumOff[cmpt][ownIndex] += abs(upper[cmpt][faceIndex])
				# sumOff[cmpt][neiIndex] += abs(lower[cmpt][faceIndex])
		
		# Ensure the matrix is diagonally dominant and that the diagonal coefficient
		# is positive
		# for cmpt in range(noComponents):
			
			# for cellIndex in range(nCells):
				
				# diag[cmpt][cellIndex] = max(abs(diag[cmpt][cellIndex]), sumOff[cmpt][cellIndex])
		
		# ... then relax
		diag /= alpha
		
		# Finally add the relaxation contribution to the source
		source += (diag - diagOld) * self.data.psi_.data.cellValues_
	
	def deRelax(self, fvSolutionParameters):
		
		alpha = 1.0
		
		try:
			alpha = fvSolutionParameters['impUR']
		except:
			None
		
		self.data.diag_ *= alpha
	
	
	def residual(self):
		
		noComponents = self.data.psi_.noComponents_
		
		residual = np.empty(noComponents, dtype = float)
		normFactor = np.full(noComponents, None)
		
		for cmpt in range(noComponents):
			
			residual[cmpt], normFactor[cmpt] =								   \
				fvSolution.calculateResidual(self, cmpt, normFactor[cmpt])
		
		return residual
	
	
	def defineRowStart(mesh):
		
		# Total number of cells
		nCells = mesh.read("geometryData.meshSize")
		# Number of internal faces
		nInternalFaces = np.size(mesh.read("connectivityData.neighbour"))
		
		owner = mesh.connectivityData.owner_
		
		# This algorithm is copied from OpenFOAM
		rowStart = np.full(nCells + 1, nInternalFaces, dtype = int)
		
		rowStart[0] = 0
		nOwnStart = 0
		i = 1
		
		for faceIndex in range(nInternalFaces):
			curOwn = owner[faceIndex]

			if (curOwn > nOwnStart):
				while (i <= curOwn):
					rowStart[i] = faceIndex
					i += 1

				nOwnStart = curOwn
		
		return rowStart
	
	
	"""------------------------ General functions ---------------------------"""
	
	def A(self):
		# Returns volScalarField object with zero gradient boundary condition -
		# first-order extrapolation to the boundary
		
		mesh = self.data.psi_.data.mesh_

		psiName = self.data.psi_.data.fieldName_
		AName = "A" + psiName
		
		# Define zero gradient boundary definition for all of the patches for diag
		psiBoundaryDefinition = self.data.psi_.data.boundaryDefinition_
		ABoundaryDefinition = dict.fromkeys(psiBoundaryDefinition)
		
		for psiPatch in psiBoundaryDefinition:
			 
			if (psiBoundaryDefinition[psiPatch][0] == 'empty'):
				ABoundaryDefinition[psiPatch] = psiBoundaryDefinition[psiPatch]
			else:
				ABoundaryDefinition[psiPatch] = 							   \
					('fixedGradient', np.zeros(1, dtype=float))
		
		AInternalField = np.zeros(1, dtype=float)
		
		AClass = volScalarField
			
		AVolField = AClass.initialize(AName, mesh, ABoundaryDefinition, AInternalField)
		
		# Copy diag reciprocal diag values to A field
		AVolField.data.cellValues_[0] = np.copy(self.data.diag_[0])
		# Convert to [s] by multiplying with volume of cells
		AVolField.data.cellValues_[0] *= mesh.geometryData.V_
		
		AVolField.setBoundaryValues(ABoundaryDefinition)
		
		return AVolField
	
	
	# Function to return the reciprocal diagonal part of the A matrix as a
	# volField, 1/A
	def rA(self):
		# Returns volScalarField object with zero gradient boundary condition -
		# first-order extrapolation to the boundary
		
		mesh = self.data.psi_.data.mesh_

		psiName = self.data.psi_.data.fieldName_
		rAName = "rA" + psiName
		
		# Define zero gradient boundary definition for all of the patches for diag
		psiBoundaryDefinition = self.data.psi_.data.boundaryDefinition_
		rABoundaryDefinition = dict.fromkeys(psiBoundaryDefinition)
		
		for psiPatch in psiBoundaryDefinition:
			 
			if (psiBoundaryDefinition[psiPatch][0] == 'empty'):
				rABoundaryDefinition[psiPatch] = psiBoundaryDefinition[psiPatch]
			else:
				rABoundaryDefinition[psiPatch] = 							   \
					('fixedGradient', np.zeros(1, dtype=float))
		
		rAInternalField = np.zeros(1, dtype=float)
		
		rAClass = volScalarField
			
		rAVolField = rAClass.initialize(rAName, mesh, rABoundaryDefinition, rAInternalField)
		
		# Copy diag reciprocal diag values to rA field
		rAVolField.data.cellValues_[0] = np.copy(np.reciprocal(self.data.diag_[0]))
		# Convert to [s] by multiplying with volume of cells
		rAVolField.data.cellValues_[0] *= mesh.geometryData.V_
		
		rAVolField.setBoundaryValues(rABoundaryDefinition)
		
		return rAVolField
	
	
	# Function to return the H() part of the matrix
	def H(self):
		# H = b - sum_N a_n * psi_n
		# Zero gradient boundary condition - first-order extrapolation to the
		# boundary
		
		mesh = self.data.psi_.data.mesh_
		
		noComponents = self.data.psi_.noComponents_
		
		# Define boundary definition for all of the patches for HField
		psiBoundaryDefinition = self.data.psi_.data.boundaryDefinition_
		HBoundaryDefinition = dict.fromkeys(psiBoundaryDefinition)
		
		for psiPatch in psiBoundaryDefinition:
			 
			if (psiBoundaryDefinition[psiPatch][0] == 'empty'):
				HBoundaryDefinition[psiPatch] = psiBoundaryDefinition[psiPatch]
			else:# (psiBoundaryDefinition[psiPatch][0] == 'fixedValue'):
				HBoundaryDefinition[psiPatch] = 							   \
					('fixedGradient', np.zeros(noComponents, dtype=float))
			# else:
				# raise RuntimeError("H for velocity boundary condition " +	   \
					# psiPatch + " has not yet been implemented!")
		
		HInternalField = np.zeros(noComponents, dtype=float)
		
		# Choose the appropriate resulting class, depending on the class of the field
		HClass = eval(self.data.psi_.className_)
		
		HVolField = HClass.initialize("H", mesh, HBoundaryDefinition, HInternalField)
		
		# Set the H values to b
		HCellValues = np.copy(self.data.source_)
		
		# Subtract a_n * psi_n from H values
		psi = self.data.psi_.data.cellValues_
		
		source = self.data.source_
		
		lower = self.data.lower_
		upper = self.data.upper_

		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		# For all components
		for cmpt in range(noComponents):
			
			# For all faces
			for faceIndex in range(neighbour.size):
				
				ownIndex = owner[faceIndex]
				neiIndex = neighbour[faceIndex]
				
				HCellValues[cmpt][ownIndex] -=								   \
					upper[cmpt][faceIndex] * psi[cmpt][neiIndex]
				
				HCellValues[cmpt][neiIndex] -=								   \
					lower[cmpt][faceIndex] * psi[cmpt][ownIndex]
		
		# Convert to [m/s2] by dividing by volume of cells
		for cmpt in range(noComponents):
			HCellValues[cmpt] /= mesh.geometryData.V_
		
		HVolField.data.cellValues_ = HCellValues
		HVolField.setBoundaryValues(HBoundaryDefinition)
		
		return HVolField
	
	
	# Prints the contributions and values for a specified cell and component
	def cellAnalysis(self, cellNo, cmpt):
		
		mesh = self.data.mesh_
		
		psi = self.data.psi_.data.cellValues_[cmpt]
		
		nInternalFaces = np.size(mesh.connectivityData.neighbour_)
		
		owner = mesh.connectivityData.owner_
		neighbour = mesh.connectivityData.neighbour_
		
		source = self.data.source_[cmpt]
		
		lower = self.data.lower_[cmpt]
		diag = self.data.diag_[cmpt]
		upper = self.data.upper_[cmpt]
		
		resultArray = np.empty(0, dtype = object)
		
		counter = 0
		
		for faceIndex in range(nInternalFaces):
			
			ownIndex = owner[faceIndex]
			neiIndex = neighbour[faceIndex]
			
			if (ownIndex == cellNo):
				
				resultArray = np.append(resultArray, None)
				
				resultArray[counter] = np.array([upper[faceIndex], neiIndex])
				
				counter += 1
			
			if (neiIndex == cellNo):
				
				resultArray = np.append(resultArray, None)
				resultArray[counter] = np.array([lower[faceIndex], ownIndex])
				
				counter += 1
		
		# Printing out results
		print()
		
		for nei in range(np.size(resultArray)):
			
			print("an to cell " + str(resultArray[nei][1]) + " is " +		   \
				str(resultArray[nei][0]) + ". The value in the cell is " + str(psi[int(resultArray[nei][1] + 1e-10)]))
		
		print("ap in cell " + str(cellNo) + " is " + str(diag[cellNo]) +	   \
			". The values in the cell is " + str(psi[cellNo]))
		
		print("b in cell " + str(cellNo) + " is " + str(source[cellNo]))
		
		print()

	"""------------------------ Defining operators --------------------------"""
	def __add__(self, other): # + operator
		
		if (self.data.mesh_ != other.data.mesh_):
			raise RuntimeError(												   \
				"The matrices are not assembled using the same mesh!")
		
		if (self.data.psi_ != other.data.psi_):
			raise RuntimeError(												   \
				"The matrices are not assembled using the same field!")
		
		# Copy
		mesh = self.data.mesh_
		
		psi = self.data.psi_
		
		rowStart = self.data.rowStart_
		
		# Add up the source and A matrix contributions
		source = np.copy(self.data.source_) + np.copy(other.data.source_)
		
		lower = np.copy(self.data.lower_) + np.copy(other.data.lower_)
		diag = np.copy(self.data.diag_) + np.copy(other.data.diag_)
		upper = np.copy(self.data.upper_) + np.copy(other.data.upper_)

		return fvMatrix(psi, source, lower, diag, upper, rowStart)
		
	def __sub__(self, other): # - operator
		
		if (self.data.mesh_ != other.data.mesh_):
			raise RuntimeError(												   \
				"The matrices are not assembled using the same mesh!")
		
		if (self.data.psi_ != other.data.psi_):
			raise RuntimeError(												   \
				"The matrices are not assembled using the same field!")
				
		# Copy
		mesh = self.data.mesh_
		
		psi = self.data.psi_
		
		rowStart = self.data.rowStart_
		
		# Subtract the source and A matrix contributions
		source = np.copy(self.data.source_) - np.copy(other.data.source_)
		
		lower = np.copy(self.data.lower_) - np.copy(other.data.lower_)
		diag = np.copy(self.data.diag_) - np.copy(other.data.diag_)
		upper = np.copy(self.data.upper_) - np.copy(other.data.upper_)
		
		return fvMatrix(psi, source, lower, diag, upper, rowStart)
	
	def __neg__(self): # - operator in front of a matrix
						
		# Copy
		mesh = self.data.mesh_
		
		psi = self.data.psi_
		
		rowStart = self.data.rowStart_
		
		# Subtract the source and A matrix contributions
		source = - np.copy(self.data.source_)
		
		lower = - np.copy(self.data.lower_)
		diag = - np.copy(self.data.diag_)
		upper = - np.copy(self.data.upper_)
		
		return fvMatrix(psi, source, lower, diag, upper, rowStart)
	
	def __eq__(self, other): # == operator
		
		if (self.data.mesh_ != other.data.mesh_):
			raise RuntimeError(												   \
				"The matrices are not assembled using the same mesh!")
		
		if (self.data.psi_ != other.data.psi_):
			raise RuntimeError(												   \
				"The matrices are not assembled using the same field!")
				
		# Copy
		mesh = self.data.mesh_
		
		psi = self.data.psi_
		
		rowStart = self.data.rowStart_
		
		# Subtract the source and A matrix contributions
		source = np.copy(self.data.source_) + np.copy(other.data.source_)
		
		lower = np.copy(self.data.lower_) - np.copy(other.data.lower_)
		diag = np.copy(self.data.diag_) - np.copy(other.data.diag_)
		upper = np.copy(self.data.upper_) - np.copy(other.data.upper_)
		
		return fvMatrix(psi, source, lower, diag, upper, rowStart)

"""
// ************************************************************************* //
"""
