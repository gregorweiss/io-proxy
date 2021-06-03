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
    Mesh data needed to do the Finite Volume discretisation.
    
    Note: Cell and face centres are calculated using averages of vertices.
    Note: The mesh can be initialized using the functions present here, or
    by using the pre-existing Ofpp package. This is modifyed by
    commenting/uncommenting the appropriate lines in the
    constructFromPolyMeshFolder function, and the import statement at the
    begining of this file
    
\*---------------------------------------------------------------------------*/
"""

import numpy as np
import subprocess
import re
# import Ofpp

class fvMesh:
	
	class fvMeshGeometryData:
		def __init__(self):
			V_ 			= None # Cell volumes, volume scalar field
			Sf_ 		= None # Face area vectors, surface vector field
			magSf_ 		= None # Mag face area vectors, surface scalar field
			C_ 			= None # Cell centres, volume vector field
			Cf_ 		= None # Face centres, surface vector field
			meshSize_	= None # Number of cells in the mesh
		
	class fvMeshConnectivityData:
		def __init__(self):
			points_ 	= None # Points list
			faces_		= None # Faces list
			owner_		= None # Owners list
			neighbour_	= None # Neighbours list
			boundary_	= None # Boundaries list
	
	"""------------------------- Constructors -------------------------------"""		
	# Main constructor
	def __init__(self, volumes, surfaces, cellCentres, faceCentres, points,    \
				 faces, owner, neighbour, boundary):
		
		self.geometryData = self.fvMeshGeometryData()
		self.connectivityData = self.fvMeshConnectivityData()
		
		# Geometry data
		self.geometryData.V_ 			= volumes
		self.geometryData.Sf_ 			= surfaces
		self.calculateMagSf()
		self.geometryData.C_ 			= cellCentres
		self.geometryData.Cf_ 			= faceCentres
		self.geometryData.meshSize_ 	= self.geometryData.V_.size
	
		# Connectivity data
		self.connectivityData.points_	= points
		self.connectivityData.faces_	= faces
		self.connectivityData.owner_	= owner
		self.connectivityData.neighbour_= neighbour
		self.connectivityData.boundary_	= boundary
		
		print("Mesh created.\n")
	
	# Constructor which reads in OpenFOAM-format polyMesh folder and calculates
	# the necessary metrics from that
	@classmethod
	def constructFromPolyMeshFolder(cls, polyMeshLocation):
		
		print("Creating mesh.\n")
		
		# """ --------- Ofpp method --------- """
		# # Reference: https://github.com/xu-xianghua/ofpp
		# mesh = Ofpp.FoamMesh(".")
		# points 			= mesh.points
		# faces 			= np.array(mesh.faces, dtype=object)
		# owner			= np.array(mesh.owner)
		# neighbour		= np.array(mesh.neighbour[:mesh.num_inner_face])
			
		# boundaryArray = np.empty(0, dtype = set)
		
		# nBoundary = 0
		# for i in mesh.boundary:
			# boundaryName = i.decode('ASCII')
			# boundaryType = (mesh.boundary[i])[0].decode('ASCII')
			# nFaces = (mesh.boundary[i])[1]
			# startFace = (mesh.boundary[i])[2]
			# boundaryArray = np.append(boundaryArray, None)
			# boundaryArray[nBoundary] =										   \
				# [boundaryName, boundaryType, nFaces, startFace]
				
			# nBoundary += 1
		
		# boundary = boundaryArray
		# """ --------- ----------- --------- """	

		""" --------- Custom method --------- """
		# Connectivity data
		points			= cls.readPointsFromPolyMesh(polyMeshLocation)
		faces 			= cls.readFacesFromPolyMesh(polyMeshLocation)
		owner			= cls.readOwnerFromPolyMesh(polyMeshLocation)
		neighbour		= cls.readNeighbourFromPolyMesh(polyMeshLocation)
		boundary 		= cls.readBoundaryFromPolyMesh(polyMeshLocation)
		""" --------- ----------- --------- """
		
		# Geometry data
		faceCentres 	= cls.calculateFaceCentres(points, faces)
		surfaces 		= cls.calculateFaceAreaVectors(points, faces,		   \
						  faceCentres)
		cellCentres 	= cls.calculateCellCentres(faceCentres, owner, 		   \
						  neighbour)
		volumes 		= cls.calculateCellVolumes(owner, neighbour, surfaces, \
						  faceCentres, cellCentres)
		
		return cls(volumes, surfaces, cellCentres, faceCentres, points, faces, \
				   owner, neighbour, boundary)
	
	"""---------------------- polyMesh reading functions --------------------"""
	
	@classmethod
	def readPointsFromPolyMesh(cls, polyMeshLocation):
		
		# Read in file
		file = open(polyMeshLocation + "/points").read()
		# Remove comments
		file = cls.removeCppComments(file)
		# Remove FoamFile dictionary
		file = cls.removeFoamFile(file)
		# Split the file into tokens
		file = file.split()
		
		nPoints = int(file[0])
		
		# Initialise the end result: array of points
		pointsArray = np.empty((nPoints, 3), dtype = float)
		
		# Variable used to know in which point the file is
		pointCounter = 0
		# Variable used to determine the X, Y or Z coordinate
		componentCounter = 0
		for i in file[2:-1]:
			coordinate = (i.replace('(', '')).replace(')', '')
			
			pointsArray[pointCounter][componentCounter] = coordinate
			
			if (componentCounter == 2):
				componentCounter = 0
				pointCounter += 1
			else:
				componentCounter += 1
				
		return pointsArray

	
	@classmethod
	def readFacesFromPolyMesh(cls, polyMeshLocation):
		
		# Read in file
		file = open(polyMeshLocation + "/faces").read()
		# Remove comments
		file = cls.removeCppComments(file)
		# Remove FoamFile dictionary
		file = cls.removeFoamFile(file)
		# Split the file into tokens
		file = file.split()
		
		nFaces = int(file[0])

		facesArray = np.empty(nFaces, dtype = object)
		
		# Variable used to know in which face the file is
		faceCounter = 0
		# Variable used to know which is the number of the point in the face
		pointInFaceCounter = 0
		for i in file[2:-1]:
			
			if "(" in i:
				splitI = i.split("(")
				# Find out the size of the new face
				newFaceSize = int(splitI[0])
				# Initialize the new face
				newFace = np.empty(newFaceSize, dtype = int)
				# Insert the first point of the new face
				newFace[0] = int(splitI[1])
				pointInFaceCounter = 0
			elif ")" in i:
				pointInFaceCounter += 1
				newFace[pointInFaceCounter] = int(i.replace(")", ""))
				facesArray[faceCounter] = newFace
				faceCounter += 1
			else:
				pointInFaceCounter += 1
				newFace[pointInFaceCounter] = int(i)

		return facesArray
	
	@classmethod
	def readOwnerFromPolyMesh(cls, polyMeshLocation):
		
		# Read in file
		file = open(polyMeshLocation + "/owner").read()
		# Remove comments
		file = cls.removeCppComments(file)
		# Remove FoamFile dictionary
		file = cls.removeFoamFile(file)
		# Split the file into tokens
		file = file.split()
		
		nOwners = int(file[0])
		
		ownerArray = np.empty(nOwners, dtype = int)
		
		ownerCounter = 0
		for i in file[2:-1]:
			ownerArray[ownerCounter] = i
			ownerCounter += 1
		
		return ownerArray
	
	@classmethod
	def readNeighbourFromPolyMesh(cls, polyMeshLocation):
		
		# Read in file
		file = open(polyMeshLocation + "/neighbour").read()
		# Remove comments
		file = cls.removeCppComments(file)
		# Remove FoamFile dictionary
		file = cls.removeFoamFile(file)
		# Split the file into tokens
		file = file.split()
		
		nNeighbours = int(file[0])
		
		neighbourArray = np.empty(nNeighbours, dtype = int)
		
		neighbourCounter = 0
		for i in file[2:-1]:
			neighbourArray[neighbourCounter] = i
			neighbourCounter += 1
		
		return neighbourArray
	
	@classmethod
	def readBoundaryFromPolyMesh(cls, polyMeshLocation):
		
		# Read in file
		file = open(polyMeshLocation + "/boundary").read()
		# Remove comments
		file = cls.removeCppComments(file)
		# Remove FoamFile dictionary
		file = cls.removeFoamFile(file)
		# Split the file into tokens
		file = file.split()
		
		nBoundary = int(file[0])

		# Initialise the end result: array of boundaries
		boundaryArray = np.empty(nBoundary, dtype = set)

		# Loop to put in place holder values into the final array
		for i in range(nBoundary):
			boundaryArray[i] = ["boundaryName", "boundaryType", 0, 0]

		# Variables/flags used for going through the boundary file
		boundaryCounter = 0
		lookingForName = True
		
		lookingForType = False
		lookingFornFaces = False
		foundnFaces = False
		lookingForstartFace = False
		foundstartFace = False
		assembledBoundary = False
		for i in file[2:-1]:
			
			if (i == "}") or (i == "{"):
				continue
			elif(lookingForName):
				boundaryName = i
				lookingForName = False
				lookingForType = True
			elif (lookingForType):
				if (i == 'type'):
					continue
				else:
					boundaryType = i.replace(";", "")
					lookingForType = False
					lookingFornFaces = True
			elif (lookingFornFaces):
				if (i == 'nFaces'):
					foundnFaces = True
				elif (not foundnFaces):
					continue
				else:
					nFaces = i.replace(";", "")
					foundnFaces = False
					lookingFornFaces = False
					lookingForstartFace = True
			elif (lookingForstartFace):
				if (i == 'startFace'):
					foundstartFace = True
				elif (not foundstartFace):
					continue
				else:
					startFace = i.replace(";", "")
					foundstartFace = False
					lookingForstartFace = False
					lookingForName = True
					
					assembledBoundary = True
			
			if (assembledBoundary):
				boundaryArray[boundaryCounter][0] = boundaryName
				boundaryArray[boundaryCounter][1] = boundaryType
				boundaryArray[boundaryCounter][2] = int(nFaces)
				boundaryArray[boundaryCounter][3] = int(startFace)
				
				boundaryCounter += 1
				assembledBoundary = False
		
		return boundaryArray
		
	def removeCppComments(text):
		def replacer(match):
			s = match.group(0)
			if s.startswith('/'):
				return " " # note: a space and not an empty string
			else:
				return s
		pattern = re.compile(
			r'//.*?$|/\*.*?\*/|\'(?:\\.|[^\\\'])*\'|"(?:\\.|[^\\"])*"',		   \
			re.DOTALL | re.MULTILINE
		)
		
		return re.sub(pattern, replacer, text)
		
	def removeFoamFile(file):
		inFoamFile = False
		lineCounter = 0
		for i in file:
			lineCounter += 1
			if i == '}':
				break

		return file[lineCounter:]
		
					
	"""------------------------- Geometric functions ------------------------"""
	
	# Function to calculate the face area magnitudes using the face area vectors
	def calculateMagSf(self):
		
		# Define the size of the result array
		self.geometryData.magSf_ = np.empty(np.size(self.geometryData.Sf_,	   \
								   axis = 0))
		
		# The result is the magnitude of face area vectors
		for i in range(self.geometryData.magSf_.size): # linalg.norm does
													   # not seem to
													   # vectorise
			self.geometryData.magSf_[i] = np.linalg.norm( 					   \
										  self.geometryData.Sf_[i])
		
	def calculateFaceAreaVectors(pointList, faceList, faceCentresList):
		
		result = np.zeros((np.size(faceList, axis = 0), 3))
		
		# For all faces
		for faceIndex in range(np.size(result, axis = 0)):
			
			# The centre of the face
			faceCentre = faceCentresList[faceIndex]
			
			# For all points in the face
			for j in range(np.size(faceList[faceIndex])):
				
				# First point:
				p1 = pointList[(faceList[faceIndex])[j]]
				
				# Second point:
				# If not at the last point in the face
				if (j < np.size(faceList[faceIndex]) - 1):
					
					p2 = pointList[(faceList[faceIndex])[j+1]]
				
				# If at the last point of the face, the next point is the first
				# (0-th)
				else:
					p2 = pointList[(faceList[faceIndex])[0]]
				
				# The vectors which define the triangle
				vec1 = p2 - p1
				vec2 = faceCentre - p1
				
				# Face area vector of the triangle between the current and
				# next points and the face centre
				result[faceIndex] += 0.5 * np.cross(vec1, vec2)
		
		return result
	
	def calculateFaceCentres(pointList, faceList):
		# Average of coordinates of vertices
		
		result = np.zeros((np.size(faceList, axis = 0), 3))
		
		# For all faces
		for i in range(np.size(faceList, axis = 0)):
			nPointsInFace = 0
			
			# For every point in face
			for j in faceList[i]:
				
				# Add up the coordinates of points
				result[i] += pointList[j]
				
				# Counter of numbers of points in face
				nPointsInFace +=1
			
			# Divide sums of coordinates by number of points in face
			result[i] /= nPointsInFace
			
		return result
		
	def calculateCellCentres(faceCentresList, ownerList, neighbourList):
		# Average of coordinates of face centres
		
		result = np.zeros((ownerList.max() + 1, 3))
		
		# Counter of number of faces in cell
		nFacesInCells = np.zeros(ownerList.max() + 1)
		
		# For all owners
		for i in range(ownerList.size):
			result[ownerList[i]] += faceCentresList[i]
			nFacesInCells[ownerList[i]] += 1
		
		# For all neighbours
		for i in range(neighbourList.size):
			result[neighbourList[i]] += faceCentresList[i]
			nFacesInCells[neighbourList[i]] += 1
		
		# Divide sums of coordinates by number of faces in cells
		for i in range(nFacesInCells.size):
			result[i] /= nFacesInCells[i]
		
		return result
		
	def calculateCellVolumes(ownerList, neighbourList, surfaces, faceCentres,  \
							 cellCentres):
		# V = sum_faces{1 / 3 * (faceCentre - cellCentre) * faceAreaVector}
		
		result = np.zeros(ownerList.max() + 1)
		
		facesInCellsList = np.empty(ownerList.max() + 1, dtype=object)
		
		# Go through owners and add faces to cells
		for i in range(ownerList.size):
			facesInCellsList[ownerList[i]] = np.append(facesInCellsList[ 	   \
			                                           ownerList[i]], [i])
			
		# Go through neighbours and add faces to cells
		for i in range(neighbourList.size):
			facesInCellsList[neighbourList[i]] = np.append(facesInCellsList[   \
												 neighbourList[i]], [i])
			
		# Delete the 0-th element in the arrays (all cell face lists start with
		# "None" without this)
		for i in range(np.size(facesInCellsList, axis = 0)):
			facesInCellsList[i] = np.delete(facesInCellsList[i], 0)
		
		# Final result calculation
		for i in range(result.size):
			for j in facesInCellsList[i]:
				result[i] += np.linalg.norm((faceCentres[j] - cellCentres[i])  \
							 * surfaces[j])
			result[i] = result[i] / 3
		
		return result

	
	"""------------------------- Access functions ---------------------------"""
	
	# Generic data access function. Input is a string which is the name of the
	# data we are interested in, et.g. 'geometryData.V',
	# 'connectivityData.owner', etc. Additionaly, the second argument can also
	# be the index of the point/face/cell, etc. we are interested in
	def read(self, inputString, i = None):
		
		readListString = "self." + inputString + "_"
		
		if (i == None):
			return eval(readListString)
		else:
			try:
				return eval(readListString + "[" + str(i) + "]")
			except RuntimeError:
				print("The requested index is out of bounds!")
				
	# Function which prints all of the data present in the mesh
	def printAll(self):		
		
		print("")
		print("List of volumes:")
		print(self.read("geometryData.V"))
		print("")
		print("List of face area vectors:")
		print(self.read("geometryData.Sf"))
		print("")
		print("List of face area vector magnitudes:")
		print(self.read("geometryData.magSf"))
		print("")
		print("List of cell centres:")
		print(self.read("geometryData.C"))
		print("")
		print("List of face centres:")
		print(self.read("geometryData.Cf"))
		print("")
		print("List of points:")
		print(self.read("connectivityData.points"))
		print("")
		print("List of faces:")
		print(self.read("connectivityData.faces"))
		print("")
		print("List of owners:")
		print(self.read("connectivityData.owner"))
		print("")
		print("List of neighbours:")
		print(self.read("connectivityData.neighbour"))
		print("")
		print("List of boundaries:")
		print(self.read("connectivityData.boundary"))
		print("")
		
		print("Number of volumes: " + 										   \
			  str(np.size(self.read("geometryData.V"), axis = 0)))
		print("Number of face area vectors: " + 							   \
			  str(np.size(self.read("geometryData.Sf"), axis = 0)))
		print("Number of face area vector magnitudes: " + 					   \
			  str(np.size(self.read("geometryData.magSf"), axis = 0)))
		print("Number of cell centres: " + 									   \
			  str(np.size(self.read("geometryData.C"), axis = 0)))
		print("Number of faceCentres: " + 									   \
			  str(np.size(self.read("geometryData.Cf"), axis = 0)))
		print("Number of points: " + 										   \
			  str(np.size(self.read("connectivityData.points"), axis = 0)))
		print("Number of faces: " + 										   \
			  str(np.size(self.read("connectivityData.faces"), axis = 0)))
		print("Number of owners: " + 										   \
			  str(np.size(self.read("connectivityData.owner"))))
		print("Number of neighbours: " + 									   \
			  str(np.size(self.read("connectivityData.neighbour"))))
		print("Number of boundaries: " + 									   \
			  str(np.size(self.read("connectivityData.boundary"), axis = 0)))

"""
// ************************************************************************* //
"""
